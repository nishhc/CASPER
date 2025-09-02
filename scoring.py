#!/usr/bin/env python3
"""
Score primer/guide pairs from a generated features CSV.

- Converts raw metrics -> per-feature 0..1 scores using base rules.
- Learns feature weights with entropy or PCA (variance / information).
- Composite score = weighted sum of feature scores (0..1). No geometric mean.
- Writes per-feature scores, weights, contributions, and the composite.

Author: you + helper
"""

from __future__ import annotations
import argparse, csv, math, sys, re, json
from typing import List, Dict, Tuple, Optional

# numpy is optional; if missing we fall back to uniform weights across kept features
try:
    import numpy as np
except Exception:
    np = None

# -------------------------
# CLI & small helpers
# -------------------------

def parse_range_2floats(text: str, name: str, lo_ok=0.0, hi_ok=1e9) -> Tuple[float, float]:
    m = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*-\s*([0-9]*\.?[0-9]+)\s*$", str(text or ""))
    if not m:
        raise argparse.ArgumentTypeError(f"{name} must be like '35-55'")
    a, b = float(m.group(1)), float(m.group(2))
    if a < lo_ok or b > hi_ok or a > b:
        raise argparse.ArgumentTypeError(f"{name} invalid: {a}-{b}")
    return (a, b)

def ffloat(d: Dict[str,str], key: str, default: float=0.0) -> float:
    v = d.get(key, "")
    if v is None: return default
    v = str(v).strip()
    if v == "" or v.upper() in {"NA","NAN","NULL"}: return default
    try:
        return float(v)
    except Exception:
        return default

def fint(d: Dict[str,str], key: str, default: int=0) -> int:
    try:
        return int(round(ffloat(d, key, float(default))))
    except Exception:
        return default

def cap01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))

def tri(x: float, lo: float, mid: float, hi: float) -> float:
    """Triangular preference: best at mid; 0 at/beyond lo and hi."""
    if x <= lo or x >= hi: return 0.0
    return (x - lo) / (mid - lo) if x < mid else (hi - x) / (hi - mid)

def inv_small(x: float, soft_ref: float) -> float:
    """Smaller is better; ≈1 when x << soft_ref, tends to 0 as x >> soft_ref."""
    return cap01(1.0 - (x / max(1e-9, soft_ref)))

def inv_pct(x: float, soft_ref: float) -> float:
    """Percent-like quantity where smaller than soft_ref is good."""
    return cap01(1.0 - (x / max(1e-9, soft_ref)))

def inv_pct_100(x: float) -> float:
    """Percent [0..100]; invert to 'goodness'."""
    return cap01(1.0 - x / 100.0)

def norm_offtarget(score: float) -> float:
    """BLAST off-target similarity in [0,1]; higher is worse -> invert."""
    return cap01(1.0 - score)

def mfe_soft_good(mfe: Optional[float]) -> float:
    """
    More negative MFE = more stable hairpin (worse).
    Map MFE in [-12 .. -2] to [0 .. 1]; >= -2 => 1.0; <= -12 => 0.0
    """
    if mfe is None: return 0.5
    if mfe >= -2.0: return 1.0
    if mfe <= -12.0: return 0.0
    return (mfe + 12.0) / 10.0

# -------------------------
# Parameters for scoring rules
# -------------------------

class Params:
    def __init__(
        self,
        gc_primer: Tuple[float,float],
        amplicon_len: Tuple[float,float],
        guide_margin_bp: float,
        tm_delta_warn: float,
        cross_dimer_thresh: float,
        self_dimer_thresh: float,
        seed_gc_cas9: Tuple[float,float],
        seed_gc_other: Tuple[float,float],
    ):
        self.gc_primer = gc_primer
        self.amplicon_len = amplicon_len
        self.guide_margin_bp = guide_margin_bp
        self.tm_delta_warn = tm_delta_warn
        self.cross_dimer_thresh = cross_dimer_thresh
        self.self_dimer_thresh = self_dimer_thresh
        self.seed_gc_cas9 = seed_gc_cas9
        self.seed_gc_other = seed_gc_other

# -------------------------
# Per-row feature scoring (0..1)
# -------------------------

def build_feature_scores_from_row(row: Dict[str,str], P: Params) -> Tuple[List[str], List[float]]:
    """
    Convert one CSV row to a consistent vector of feature scores in [0,1].
    Returns (feature_names, scores).
    All features here are *goodness* (higher = better).
    """
    names: List[str] = []
    vals:  List[float] = []

    # --- Specificity (invert BLAST similarity) ---
    names += ["spec_fp", "spec_rp", "spec_crrna"]
    vals  += [
        norm_offtarget(ffloat(row, "fp_offtarget_score", 0.0)),
        norm_offtarget(ffloat(row, "rp_offtarget_score", 0.0)),
        norm_offtarget(ffloat(row, "crrna_offtarget_score", 0.0)),
    ]

    # Penalize exact full-length off-target matches (1=has exact -> bad)
    # turn into a goodness (no exact=1, exact=0)
    names += ["spec_no_exact_fp", "spec_no_exact_rp", "spec_no_exact_crrna"]
    vals  += [
        1.0 - cap01(ffloat(row, "fp_has_exact_full_hit", 0.0)),
        1.0 - cap01(ffloat(row, "rp_has_exact_full_hit", 0.0)),
        1.0 - cap01(ffloat(row, "crrna_has_exact_full_hit", 0.0)),
    ]

    # --- Dimerization / ΔTm ---
    names += ["cross_dimer", "fp_self_dimer", "rp_self_dimer", "delta_tm"]
    vals  += [
        inv_pct(ffloat(row, "fp_rp_cross_dimer_pct", 0.0), P.cross_dimer_thresh),
        inv_pct(ffloat(row, "fp_self_dimer_pct", 0.0),     P.self_dimer_thresh),
        inv_pct(ffloat(row, "rp_self_dimer_pct", 0.0),     P.self_dimer_thresh),
        inv_small(ffloat(row, "delta_tm_C", 0.0),          P.tm_delta_warn),
    ]

    # 3′ complementarity runs; shorter is better
    names += ["fp_rp_3p_cross_run", "fp_3p_self_run", "rp_3p_self_run"]
    vals  += [
        inv_small(ffloat(row, "fp_rp_3p_cross_run", 0.0), 4.0),
        inv_small(ffloat(row, "fp_3p_self_run", 0.0),     4.0),
        inv_small(ffloat(row, "rp_3p_self_run", 0.0),     4.0),
    ]

    # 3′ G-quad flags (0 good, 1 bad -> invert)
    names += ["fp_no_gquad3p", "rp_no_gquad3p"]
    vals  += [1.0 - cap01(ffloat(row, "fp_gquad_3p", 0.0)),
              1.0 - cap01(ffloat(row, "rp_gquad_3p", 0.0))]

    # --- Geometry ---
    lm = ffloat(row, "left_margin_bp", 0.0)
    rm = ffloat(row, "right_margin_bp", 0.0)
    names += ["left_margin", "right_margin", "amplicon_len"]
    vals  += [
        min(1.0, lm / max(1.0, P.guide_margin_bp)),
        min(1.0, rm / max(1.0, P.guide_margin_bp)),
        tri(ffloat(row, "amplicon_len", 0.0), P.amplicon_len[0], (P.amplicon_len[0]+P.amplicon_len[1])/2.0, P.amplicon_len[1]),
    ]

    # On-target amplicon uniqueness (exactly 1 is ideal)
    names += ["amplicon_unique"]
    vals  += [1.0 if fint(row, "amplicon_hits_on_target", 1) == 1 else 0.0]

    # --- Composition ---
    lo_gc_p, hi_gc_p = P.gc_primer
    mid_gc_p = (lo_gc_p + hi_gc_p) / 2.0
    names += ["fp_gc", "rp_gc"]
    vals  += [
        tri(ffloat(row, "fp_gc_pct", 0.0), lo_gc_p, mid_gc_p, hi_gc_p),
        tri(ffloat(row, "rp_gc_pct", 0.0), lo_gc_p, mid_gc_p, hi_gc_p),
    ]

    # Seed GC (Cas-specific preferences)
    seed_gc = ffloat(row, "seed_gc_pct", 50.0)
    cas = (row.get("cas_type","") or "").lower()
    lo_seed, hi_seed = (P.seed_gc_cas9 if cas == "cas9" else P.seed_gc_other)
    names += ["seed_gc"]
    vals  += [tri(seed_gc, lo_seed, (lo_seed+hi_seed)/2.0, hi_seed)]

    # Seed homopolymer (0/1; want 0)
    names += ["seed_no_hpoly"]
    vals  += [1.0 - cap01(ffloat(row, "seed_has_homopolymer", 0.0))]

    # Global/primer complexity (short homopolymers are good)
    names += ["amplicon_low_hpoly", "fp_low_hpoly", "rp_low_hpoly"]
    vals  += [
        inv_small(ffloat(row, "amplicon_max_run", 1.0), 8.0),
        inv_small(ffloat(row, "fp_max_hpoly_run", 1.0), 6.0),
        inv_small(ffloat(row, "rp_max_hpoly_run", 1.0), 6.0),
    ]

    # Cross-hybridization with guide (lower % complementarity is better)
    names += ["fp_low_guide_cross", "rp_low_guide_cross"]
    vals  += [
        inv_pct_100(ffloat(row, "fp_guide_cross_pct", 0.0)),
        inv_pct_100(ffloat(row, "rp_guide_cross_pct", 0.0)),
    ]

    # Guide accessibility / secondary structure
    names += ["guide_seed_access", "guide_soft_mfe"]
    vals  += [
        cap01(ffloat(row, "guide_seed_unpaired_frac", 0.5)),
        mfe_soft_good(ffloat(row, "guide_mfe_kcal", 0.0)),
    ]

    # Ensure bounded [0,1]
    vals = [cap01(v) for v in vals]
    return names, vals

# -------------------------
# Weight learning (entropy / PCA)
# -------------------------

def entropy_weights(X: "np.ndarray") -> "np.ndarray":
    """
    Entropy weighting: columns (features) with higher dispersion/information get larger weights.
    """
    eps = 1e-12
    # Normalize each column to probability-like distribution across rows
    col_sum = X.sum(axis=0)
    col_sum[col_sum == 0] = eps
    P = X / col_sum
    n = max(2, X.shape[0])
    k = 1.0 / math.log(n)
    H = -k * (P * np.log(P + eps)).sum(axis=0)   # entropy per feature
    D = 1.0 - H                                  # divergence (information)
    D[D < 0] = 0
    s = D.sum()
    return D / (s if s > eps else 1.0)

def pca_weights(X: "np.ndarray") -> "np.ndarray":
    """
    PCA-based weights: squared loadings integrated over PCs, weighted by explained variance.
    """
    # standardize columns
    X0 = X - X.mean(axis=0, keepdims=True)
    s = X0.std(axis=0, ddof=1, keepdims=True); s[s == 0] = 1.0
    Z = X0 / s
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    expl = (S**2) / max(1.0, Z.shape[0]-1)      # eigenvalues
    expl = expl / max(1e-12, expl.sum())        # variance fractions
    load = Vt.T                                  # feature x PC
    importance = (load**2 * expl).sum(axis=1)    # sum over PCs
    importance[importance < 0] = 0
    s2 = importance.sum()
    return importance / (s2 if s2 > 1e-12 else 1.0)

def learn_feature_weights(rows: List[Dict[str,str]], P: Params, mode: str):
    # Build full score matrix
    name_ref, v0 = build_feature_scores_from_row(rows[0], P)
    mat = [v0]
    for r in rows[1:]:
        _, v = build_feature_scores_from_row(r, P)
        mat.append(v)

    if np is None:
        # Fallback: uniform weights across non-constant features
        # (and zero weight for constant features)
        # Detect variance
        m = len(mat)
        k = len(name_ref)
        means = [sum(row[j] for row in mat)/m for j in range(k)]
        vars_ = []
        for j in range(k):
            v = sum((row[j]-means[j])**2 for row in mat)/(m-1 if m>1 else 1)
            vars_.append(v)
        keep = [i for i, vv in enumerate(vars_) if vv > 1e-9]
        w = [0.0]*k
        if keep:
            unif = 1.0/len(keep)
            for i in keep:
                w[i] = unif
        return name_ref, w

    X = np.array(mat, dtype=float)  # rows x features
    # Drop (or zero-weight) features with ~zero variance
    var = X.var(axis=0)
    keep_idx = [i for i, vv in enumerate(var) if vv > 1e-9]
    if not keep_idx:
        # all constant → uniform anyway
        return name_ref, [1.0/len(name_ref)]*len(name_ref)

    Xk = X[:, keep_idx]
    if mode == "pca":
        wk = pca_weights(Xk)
    else:
        wk = entropy_weights(Xk)

    # Place back into full vector
    w_full = np.zeros(X.shape[1], dtype=float)
    w_full[keep_idx] = wk
    w_full = w_full / max(1e-12, w_full.sum())
    return name_ref, w_full.tolist()

# -------------------------
# Scoring workflow
# -------------------------

def score_rows(rows: List[Dict[str,str]], P: Params, mode: str):
    # Learn weights on *scored* features
    feat_names, weights = learn_feature_weights(rows, P, mode)
    # Apply to all rows
    out_rows = []
    for r in rows:
        names, s = build_feature_scores_from_row(r, P)
        # Safety: order must match feat_names
        if names != feat_names:
            raise RuntimeError("Feature ordering mismatch; check your CSV columns.")
        # Composite = weighted *sum* of scores (weights sum to 1) → 0..1
        composite = sum(w * si for w, si in zip(weights, s))
        # Build outputs
        r_out = dict(r)  # start with original columns
        r_out["score_mode"] = mode
        r_out["composite_score"] = f"{composite:.4f}"
        # Per-feature columns: scores, weights, contributions
        for n, si, wi in zip(names, s, weights):
            r_out[f"score_{n}"] = f"{si:.4f}"
            r_out[f"weight_{n}"] = f"{wi:.6f}"
            r_out[f"contrib_{n}"] = f"{(wi*si):.6f}"
        out_rows.append(r_out)
    return feat_names, weights, out_rows

# -------------------------
# I/O
# -------------------------

def load_csv(path: str) -> List[Dict[str,str]]:
    rows: List[Dict[str,str]] = []
    with open(path, newline="") as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            rows.append(r)
    if not rows:
        raise RuntimeError("No rows in input CSV.")
    return rows

def write_csv(path: str, rows: List[Dict[str,str]], extra_order: List[str]):
    # Build field order: original fields + our extras grouped
    # Infer original header from first row (minus our extras)
    base_keys = [k for k in rows[0].keys() if not (k.startswith("score_") or k.startswith("weight_") or k.startswith("contrib_") or k in {"score_mode","composite_score"})]
    # Compose final header
    fields = list(base_keys) + ["score_mode","composite_score"] + extra_order
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def build_arg_parser():
    p = argparse.ArgumentParser(description="Score primer/guide pairs from a features CSV (no geometric mean).")
    p.add_argument("--input-csv", required=True, help="features CSV from your generator")
    p.add_argument("--output-csv", required=True, help="output scored CSV")
    p.add_argument("--mode", choices=["entropy","pca"], default="entropy", help="weight learning mode")
    # base rules (tune if needed)
    p.add_argument("--gc-primer", default="30-70", help="primer GC%% optimal range lo-hi (default 30-70)")
    p.add_argument("--amplicon-len", default="100-250", help="amplicon length optimal range lo-hi (bp)")
    p.add_argument("--guide-margin", type=float, default=20.0, help="desired min margin from guide to each primer end (bp)")
    p.add_argument("--tm-delta", type=float, default=10.0, help="ΔTm soft threshold where larger is worse (°C)")
    p.add_argument("--cross-dimer-thresh", type=float, default=40.0, help="cross-dimer %% soft threshold")
    p.add_argument("--self-dimer-thresh", type=float, default=40.0, help="self-dimer %% soft threshold")
    p.add_argument("--seed-gc-cas9", default="40-60", help="seed GC%% optimal range for Cas9 guides (lo-hi)")
    p.add_argument("--seed-gc-other", default="35-55", help="seed GC%% optimal range for Cas12a/Cas13 (lo-hi)")
    return p

def main(argv=None):
    ap = build_arg_parser()
    args = ap.parse_args(argv)

    rows = load_csv(args.input_csv)

    P = Params(
        gc_primer=parse_range_2floats(args.gc_primer, "gc-primer"),
        amplicon_len=parse_range_2floats(args.amplicon_len, "amplicon-len"),
        guide_margin_bp=float(args.guide_margin),
        tm_delta_warn=float(args.tm_delta),
        cross_dimer_thresh=float(args.cross_dimer_thresh),
        self_dimer_thresh=float(args.self_dimer_thresh),
        seed_gc_cas9=parse_range_2floats(args.seed_gc_cas9, "seed-gc-cas9"),
        seed_gc_other=parse_range_2floats(args.seed_gc_other, "seed-gc-other"),
    )

    # Score
    feat_names, weights, out_rows = score_rows(rows, P, args.mode)

    # Build an extra ordered block for per-feature results
    extra_cols: List[str] = []
    for n in feat_names:
        extra_cols += [f"weight_{n}"]
    for n in feat_names:
        extra_cols += [f"score_{n}"]
    for n in feat_names:
        extra_cols += [f"contrib_{n}"]

    write_csv(args.output_csv, out_rows, extra_cols)

    # Short console summary
    print(f"[OK] Scored {len(out_rows)} rows with mode={args.mode}.")
    print("Top 10 features by learned weight:")
    top = sorted(zip(feat_names, weights), key=lambda z: z[1], reverse=True)[:10]
    for n, w in top:
        print(f"  {n:24s}  w={w:.4f}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)
