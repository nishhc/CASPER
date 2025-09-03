#!/usr/bin/env python3
from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Tuple, Optional, List, Dict, Union
import argparse
import os
import sys
import re
import csv
import tempfile
import subprocess
import json

# =========================
# Small utils
# =========================

def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")
    return seq.translate(table)[::-1]

def gc_pct(seq: str) -> float:
    s = seq.upper()
    g = s.count("G"); c = s.count("C")
    atgc = sum(s.count(x) for x in "ACGT")
    return 100.0 * (g + c) / max(1, atgc)

def max_homopolymer_run(seq: str) -> int:
    if not seq:
        return 0
    run = best = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            run += 1
            if run > best: best = run
        else:
            run = 1
    return best

def has_homopolymer(seq: str, max_run: int) -> bool:
    return max_homopolymer_run(seq) >= max_run

def wallace_tm(seq: str) -> float:
    s = seq.upper()
    a = s.count("A"); t = s.count("T"); g = s.count("G"); c = s.count("C")
    return 2.0*(a+t) + 4.0*(g+c)

def gentle_gc_clamp_ok(seq: str) -> int:
    tail = seq[-3:] if len(seq) >= 3 else seq
    return tail.count("G") + tail.count("C")

def find_all_indices(hay: str, needle: str) -> List[int]:
    out = []
    i = 0
    L = len(needle)
    if L == 0:
        return out
    while True:
        j = hay.find(needle, i)
        if j == -1:
            break
        out.append(j)
        i = j + 1  # overlapping allowed
    return out

def hamming(a: str, b: str) -> int:
    assert len(a) == len(b)
    return sum(1 for x,y in zip(a,b) if x != y)

def count_approx_matches(hay: str, pattern: str, max_mm: int, extra_string: Optional[str]=None) -> int:
    L = len(pattern)
    n = 0
    for i in range(0, len(hay) - L + 1):
        if hamming(hay[i:i+L], pattern) <= max_mm:
            n += 1
    if extra_string:
        p = extra_string
        for i in range(0, len(hay) - L + 1):
            if hamming(hay[i:i+L], p) <= max_mm:
                n += 1
    return n

def max_complementarity_percent(a: str, b: str) -> float:
    comp = str.maketrans("ACGT", "TGCA")
    aU = a.upper()
    bC = b.upper().translate(comp)[::-1]  # reverse-complement for pairing
    la, lb = len(aU), len(bC)
    if la == 0 or lb == 0:
        return 0.0
    max_match = 0
    for shift in range(-(lb-1), la):
        matches = 0
        overlap = 0
        for i in range(lb):
            j = i + shift
            if 0 <= j < la:
                overlap += 1
                if bC[i] == aU[j]:
                    matches += 1
        if overlap > 0:
            max_match = max(max_match, matches)
    denom = min(la, lb)
    return 100.0 * (max_match / max(1, denom))

def three_prime_run(a: str, b: str) -> int:
    aU = a.upper()
    brc = revcomp(b.upper())
    m = min(len(aU), len(brc))
    run = 0
    for k in range(1, m+1):
        if aU[-k] == brc[-k]:
            run += 1
        else:
            break
    return run

def three_prime_self_hairpin_run(seq: str) -> int:
    s = seq.upper()
    for k in range(min(12, len(s)), 0, -1):  # check up to 12 nt tail
        tail = s[-k:]
        tail_rc = revcomp(tail)
        if tail_rc in s[:-k]:
            return k
    return 0

def has_gquad_seed_3p(seq: str, window: int = 12) -> bool:
    tail = seq[-window:].upper() if len(seq) >= window else seq.upper()
    return "GGGG" in tail

def which(cmd: str) -> Optional[str]:
    from shutil import which as _which
    return _which(cmd)

# ---------- objective transforms ----------

def _tri_range(x: float, ideal: float, tol: float, hard_lo: float, hard_hi: float) -> float:
    if x <= hard_lo or x >= hard_hi:
        return 0.0
    if x == ideal:
        return 1.0
    if x < ideal:
        return max(0.0, 1.0 - (ideal - x) / max(1e-9, (ideal - hard_lo if tol == 0 else tol)))
    return max(0.0, 1.0 - (x - ideal) / max(1e-9, (hard_hi - ideal if tol == 0 else tol)))

def _band(x: float, lo: float, hi: float, soft_margin: float = 5.0) -> float:
    if lo <= x <= hi:
        return 1.0
    if x < lo:
        return max(0.0, 1.0 - (lo - x) / max(1e-9, soft_margin))
    return max(0.0, 1.0 - (x - hi) / max(1e-9, soft_margin))

def _less_is_better(x: float, good: float, bad: float) -> float:
    if x <= good: return 1.0
    if x >= bad: return 0.0
    return 1.0 - (x - good)/max(1e-9, bad - good)

def _more_is_better(x: float, good: float, cap: float) -> float:
    if x >= good: return 1.0
    return min(1.0, max(0.0, x / max(1e-9, good)))

def _mean(lst: List[float]) -> float:
    return sum(lst)/len(lst) if lst else 0.5

# =========================
# Config and I/O
# =========================

@dataclass
class Config:
    cas_type: str
    top_k: int

    primer_len_range: Tuple[int, int]
    amplicon_len_range: Tuple[int, int]
    gc_primer_range: Tuple[float, float]
    max_homopolymer: int

    crrna_len: int
    gc_crrna_range: Tuple[float, float]
    max_homopolymer_crrna: int

    cross_dimer_thresh_pct: float
    self_dimer_thresh_pct: float
    tm_delta_warn: float

    blast_word_size: int
    blast_evalue: float
    blast_task: str
    threads: int

    target_name: str
    background_fasta: Optional[str]

    guide_margin_bp: int = 20

    on_target_k: int = 12
    offtarget_max_allow: float = 0.80

def load_first_fasta_seq(path: str) -> Tuple[str, str]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA not found: {path}")
    name = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].strip() or "target"
                else:
                    break
            else:
                seq_chunks.append(line)
    if name is None:
        raise ValueError(f"No FASTA header (>) found in: {path}")
    seq = "".join(seq_chunks).upper().replace("U", "T")
    if not seq:
        raise ValueError(f"No sequence found under first header in: {path}")
    return name, seq

def parse_range_2ints(text: str, name: str) -> Tuple[int, int]:
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", text)
    if not m:
        raise argparse.ArgumentTypeError(f"{name} must be like '30-35'")
    a, b = int(m.group(1)), int(m.group(2))
    if a < 1 or b < 1 or a > b:
        raise argparse.ArgumentTypeError(f"{name} invalid: {a}-{b}")
    return (a, b)

def parse_range_2floats(text: str, name: str, lo_ok=0.0, hi_ok=100.0) -> Tuple[float, float]:
    m = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*-\s*([0-9]*\.?[0-9]+)\s*$", text)
    if not m:
        raise argparse.ArgumentTypeError(f"{name} must be like '35-55'")
    a, b = float(m.group(1)), float(m.group(2))
    if a < lo_ok or b > hi_ok or a > b:
        raise argparse.ArgumentTypeError(f"{name} invalid: {a}-{b}")
    return (a, b)

def validate_cas(t: str) -> str:
    t2 = t.strip().lower()
    if t2 not in {"cas12a", "cas9", "cas13"}:
        raise argparse.ArgumentTypeError("cas-type must be one of: cas12a, cas9, cas13")
    return t2

def default_crrna_len_for(cas_type: str) -> int:
    return {"cas12a": 23, "cas9": 20, "cas13": 28}[cas_type]

# ---------- MSA loader for conservation ----------

def load_msa_conservation(msa_path: Optional[str]) -> Optional[Dict[str, List[float]]]:
    if not msa_path:
        return None
    if not os.path.exists(msa_path):
        raise FileNotFoundError(f"MSA FASTA not found: {msa_path}")

    aln = []
    with open(msa_path) as f:
        name, seq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    aln.append("".join(seq).upper())
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if name is not None:
            aln.append("".join(seq).upper())

    if not aln or any(len(s) != len(aln[0]) for s in aln):
        raise ValueError("All MSA sequences must have identical aligned length.")

    ref = aln[0]
    ref_map = []
    ref_idx = -1
    for col, ch in enumerate(ref):
        if ch != "-":
            ref_idx += 1
            ref_map.append(col)

    cons_cols = [0.5]*len(ref)
    for col in range(len(ref)):
        counts = {"A":0,"C":0,"G":0,"T":0}
        n = 0
        for s in aln:
            b = s[col]
            if b in "ACGT":
                counts[b] += 1
                n += 1
        if n >= 3:
            cons_cols[col] = max(counts.values())/n
        else:
            cons_cols[col] = 0.5
    return {"cons_cols": cons_cols, "ref_map": ref_map}

# ---------- CLI ----------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Design FP/RP primers around CRISPR guide + features + optional background BLAST + conservation + PAM-aware scoring, with ranking."
    )
    # Required
    p.add_argument("--target-fasta", type=str, required=True)
    p.add_argument("--cas-type", type=validate_cas, required=True)
    p.add_argument("--background-fasta", type=str, default=None,
                   help="Optional background FASTA for BLAST specificity. If omitted, BLAST is skipped.")
    p.add_argument("--top-k", type=int, default=10)

    # Primer & amplicon
    p.add_argument("--primer-len", type=lambda s: parse_range_2ints(s, "primer-len"), default="30-35")
    p.add_argument("--amplicon-len", type=lambda s: parse_range_2ints(s, "amplicon-len"), default="100-250")
    p.add_argument("--gc-primer", type=lambda s: parse_range_2floats(s, "gc-primer"), default="30-70")
    p.add_argument("--max-homopolymer", type=int, default=5)

    # crRNA
    p.add_argument("--crrna-len", type=int, default=None)
    p.add_argument("--gc-crrna", type=lambda s: parse_range_2floats(s, "gc-crrna"), default="35-55")
    p.add_argument("--max-homopolymer-crrna", type=int, default=5)

    # Structure thresholds
    p.add_argument("--cross-dimer-thresh", type=float, default=40.0)
    p.add_argument("--self-dimer-thresh", type=float, default=40.0)
    p.add_argument("--tm-delta-warn", type=float, default=10.0)

    # BLAST
    p.add_argument("--blast-word-size", type=int, default=4)
    p.add_argument("--blast-evalue", type=float, default=1000.0)
    p.add_argument("--blast-task", type=str, default="blastn")
    p.add_argument("--threads", type=int, default=1)

    # Geometry
    p.add_argument("--guide-margin", type=int, default=20)

    # On-target k & specificity scale
    p.add_argument("--on-target-k", type=int, default=12)
    p.add_argument("--offtarget-max-allow", type=float, default=0.80)

    # Output
    p.add_argument("--write-csv", type=str, default="pairs_features.csv")

    # MSA & weights
    p.add_argument("--msa-fasta", type=str, default=None,
                   help="Optional MSA FASTA where the FIRST sequence is the gapped target. Used to compute conservation.")
    p.add_argument("--weights-json", type=str, default=None,
                   help="Optional JSON file with scoring weights (keys: uniq, thermo, spec, struct, geo, guide, ctx, cons, pam).")
    return p

def step1_get_config_and_sequences(argv: Optional[List[str]] = None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    target_name, target_seq = load_first_fasta_seq(args.target_fasta)

    bg_path: Optional[str] = None
    if args.background_fasta:
        if not os.path.exists(args.background_fasta):
            parser.error(f"--background-fasta not found: {args.background_fasta}")
        bg_path = os.path.abspath(args.background_fasta)

    cas = args.cas_type
    crrna_len = args.crrna_len if args.crrna_len is not None else default_crrna_len_for(cas)

    cfg = Config(
        cas_type=cas,
        top_k=max(1, int(args.top_k)),

        primer_len_range=args.primer_len,
        amplicon_len_range=args.amplicon_len,
        gc_primer_range=args.gc_primer,
        max_homopolymer=max(2, int(args.max_homopolymer)),

        crrna_len=int(crrna_len),
        gc_crrna_range=args.gc_crrna,
        max_homopolymer_crrna=max(2, int(args.max_homopolymer_crrna)),

        cross_dimer_thresh_pct=float(args.cross_dimer_thresh),
        self_dimer_thresh_pct=float(args.self_dimer_thresh),
        tm_delta_warn=float(args.tm_delta_warn),

        blast_word_size=int(args.blast_word_size),
        blast_evalue=float(args.blast_evalue),
        blast_task=str(args.blast_task),
        threads=max(1, int(args.threads)),

        target_name=target_name,
        background_fasta=bg_path,

        guide_margin_bp=max(0, int(args.guide_margin)),
        on_target_k=max(6, int(args.on_target_k)),
        offtarget_max_allow=float(args.offtarget_max_allow),
    )

    msa_info = load_msa_conservation(getattr(args, "msa_fasta", None))
    weights = None
    if getattr(args, "weights_json", None):
        with open(args.weights_json) as f:
            weights = json.load(f)

    print("\n[CONFIG]")
    for k, v in asdict(cfg).items():
        print(f"  {k}: {v}")
    print(f"\n[INPUT] target {target_name}: {len(target_seq)} bp")
    print(f"[INPUT] background FASTA: {cfg.background_fasta or '(none)'}")
    print(f"[INPUT] MSA FASTA: {getattr(args, 'msa_fasta', None) or '(none)'}")
    print(f"[INPUT] weights JSON: {getattr(args, 'weights_json', None) or '(default)'}\n")

    return cfg, target_name, target_seq, cfg.background_fasta, args.write_csv, msa_info, weights

# =========================
# Guide discovery
# =========================

@dataclass
class Guide:
    cas_type: str
    strand: str
    pam_seq: str
    pam_start: int
    protospacer_seq: str
    protospacer_start: int
    gc_pct: float
    has_homopolymer: bool

def find_crrna_sites_cas12a(target: str, guide_len: int) -> List[Guide]:
    t = target.upper()
    L = len(t)
    out: List[Guide] = []
    # + strand PAM: TTTV
    for i in range(0, L - 4 - guide_len + 1):
        pam = t[i:i+4]
        if pam.startswith("TTT") and pam[3] in "ACGT":
            prot_start = i + 4
            prot = t[prot_start: prot_start + guide_len]
            if len(prot) == guide_len:
                out.append(Guide(
                    cas_type="cas12a", strand="+", pam_seq=pam, pam_start=i,
                    protospacer_seq=prot, protospacer_start=prot_start,
                    gc_pct=gc_pct(prot), has_homopolymer=False
                ))
    # - strand PAM: AAAN on + strand
    for i in range(0, L - 4):
        pam_plus = t[i:i+4]
        if pam_plus.startswith("AAA") and pam_plus[3] in "ACGT":
            prot_end_plus = i
            prot_start_plus = prot_end_plus - guide_len
            if prot_start_plus >= 0:
                prot_plus = t[prot_start_plus: prot_end_plus]
                prot = revcomp(prot_plus)
                out.append(Guide(
                    cas_type="cas12a", strand="-", pam_seq=revcomp(pam_plus), pam_start=i,
                    protospacer_seq=prot, protospacer_start=prot_start_plus,
                    gc_pct=gc_pct(prot), has_homopolymer=False
                ))
    return out

def find_crrna_sites_cas9(target: str, guide_len: int) -> List[Guide]:
    t = target.upper()
    L = len(t)
    out: List[Guide] = []
    # + strand NGG
    for i in range(0, L - 2):
        pam = t[i:i+3]
        if len(pam) == 3 and pam[1:] == "GG":
            prot_end = i
            prot_start = prot_end - guide_len
            if prot_start >= 0:
                prot = t[prot_start: prot_end]
                out.append(Guide(
                    cas_type="cas9", strand="+", pam_seq=pam, pam_start=i,
                    protospacer_seq=prot, protospacer_start=prot_start,
                    gc_pct=gc_pct(prot), has_homopolymer=False
                ))
    # - strand CCN
    for i in range(0, L - 2):
        pam_plus = t[i:i+3]
        if len(pam_plus) == 3 and pam_plus[:2] == "CC":
            prot_start_plus = i + 3
            prot_end_plus = prot_start_plus + guide_len
            if prot_end_plus <= L:
                prot_plus = t[prot_start_plus: prot_end_plus]
                prot = revcomp(prot_plus)
                out.append(Guide(
                    cas_type="cas9", strand="-", pam_seq=revcomp(pam_plus), pam_start=i,
                    protospacer_seq=prot, protospacer_start=prot_start_plus,
                    gc_pct=gc_pct(prot), has_homopolymer=False
                ))
    return out

def find_crrna_sites_cas13(target: str, guide_len: int) -> List[Guide]:
    t = target.upper()
    L = len(t)
    out: List[Guide] = []
    for i in range(0, L - guide_len + 1):
        prot = t[i:i+guide_len]
        out.append(Guide(
            cas_type="cas13", strand="+", pam_seq="", pam_start=i,
            protospacer_seq=prot, protospacer_start=i,
            gc_pct=gc_pct(prot), has_homopolymer=False
        ))
    return out

def step2_discover_crrnas(target_seq: str, cfg: Config) -> List[Guide]:
    cas = cfg.cas_type
    if cas == "cas12a":
        raw = find_crrna_sites_cas12a(target_seq, cfg.crrna_len)
    elif cas == "cas9":
        raw = find_crrna_sites_cas9(target_seq, cfg.crrna_len)
    else:
        raw = find_crrna_sites_cas13(target_seq, cfg.crrna_len)

    lo, hi = cfg.gc_crrna_range
    filtered: List[Guide] = []
    for g in raw:
        hp = has_homopolymer(g.protospacer_seq, cfg.max_homopolymer_crrna)
        g.has_homopolymer = hp
        if hp:
            continue
        if not (lo <= g.gc_pct <= hi):
            continue
        filtered.append(g)

    print(f"[CRRNA] Found {len(raw)} raw {cas} sites; {len(filtered)} passed filters "
          f"(GC {lo}–{hi}%, homopolymer<{cfg.max_homopolymer_crrna}).")

    for g in filtered[:min(10, len(filtered))]:
        pam_info = f"PAM={g.pam_seq} at {g.pam_start}" if g.pam_seq else "PAM=NA"
        print(f"  {g.cas_type} {g.strand}  {pam_info}  prot@{g.protospacer_start}  "
              f"GC={g.gc_pct:.1f}%  seq={g.protospacer_seq}")

    return filtered

# =========================
# Primer candidates
# =========================

@dataclass
class Primer:
    start: int
    length: int
    seq: str
    strand: str           # '+' FP, '-' RP
    gc_pct: float
    has_homopolymer: bool
    tail_gc_count: int

def candidate_primers_plus_strand(target: str, cfg: Config) -> List[Primer]:
    t = target.upper()
    out: List[Primer] = []
    lo_gc, hi_gc = cfg.gc_primer_range
    for L in range(cfg.primer_len_range[0], cfg.primer_len_range[1]+1):
        for i in range(0, len(t) - L + 1):
            s = t[i:i+L]
            gc = gc_pct(s)
            if not (lo_gc <= gc <= hi_gc):
                continue
            if has_homopolymer(s, cfg.max_homopolymer):
                continue
            clamp = gentle_gc_clamp_ok(s)
            if clamp != 1:
                continue
            out.append(Primer(start=i, length=L, seq=s, strand='+',
                              gc_pct=gc, has_homopolymer=False, tail_gc_count=clamp))
    return out

def candidate_primers_minus_strand(target: str, cfg: Config) -> List[Primer]:
    t = target.upper()
    out: List[Primer] = []
    lo_gc, hi_gc = cfg.gc_primer_range
    for L in range(cfg.primer_len_range[0], cfg.primer_len_range[1]+1):
        for i in range(0, len(t) - L + 1):
            win = t[i:i+L]
            s = revcomp(win)
            gc = gc_pct(s)
            if not (lo_gc <= gc <= hi_gc):
                continue
            if has_homopolymer(s, cfg.max_homopolymer):
                continue
            clamp = gentle_gc_clamp_ok(s)
            if clamp != 1:
                continue
            out.append(Primer(start=i, length=L, seq=s, strand='-',
                              gc_pct=gc, has_homopolymer=False, tail_gc_count=clamp))
    return out

# =========================
# Pairing
# =========================

@dataclass
class PairedSet:
    guide_index: int
    guide: Guide
    fp: Primer
    rp: Primer
    amplicon_len: int

    metrics: Dict[str, float]
    flags: Dict[str, str]

def pair_primers_around_guides(target: str,
                               guides: List[Guide],
                               fps: List[Primer],
                               rps: List[Primer],
                               cfg: Config) -> List[PairedSet]:
    pairs: List[PairedSet] = []
    min_amp, max_amp = cfg.amplicon_len_range
    margin = cfg.guide_margin_bp

    for gi, g in enumerate(guides):
        left_ok = [fp for fp in fps if (fp.start + fp.length) <= (g.protospacer_start - margin)]
        guide_end = g.protospacer_start + cfg.crrna_len
        right_ok = [rp for rp in rps if rp.start >= (guide_end + margin)]

        for fp in left_ok:
            amplicon_start = fp.start
            for rp in right_ok:
                amplicon_end = rp.start + rp.length  # exclusive
                amplicon = amplicon_end - amplicon_start
                if amplicon < min_amp or amplicon > max_amp:
                    continue
                left_gap = g.protospacer_start - amplicon_start
                right_gap = amplicon_end - (g.protospacer_start + cfg.crrna_len)
                if left_gap < margin or right_gap < margin:
                    continue
                pairs.append(PairedSet(
                    guide_index=gi, guide=g, fp=fp, rp=rp,
                    amplicon_len=amplicon, metrics={}, flags={}
                ))
    return pairs

# =========================
# Feature annotations
# =========================

def annotate_structure_and_tm(ps: PairedSet, cfg: Config):
    tm_fp = wallace_tm(ps.fp.seq)
    tm_rp = wallace_tm(ps.rp.seq)
    delta_tm = abs(tm_fp - tm_rp)

    fp_self = max_complementarity_percent(ps.fp.seq, ps.fp.seq)
    rp_self = max_complementarity_percent(ps.rp.seq, ps.rp.seq)
    fp_rp_cross = max_complementarity_percent(ps.fp.seq, ps.rp.seq)

    fp_3p_self = three_prime_self_hairpin_run(ps.fp.seq)
    rp_3p_self = three_prime_self_hairpin_run(ps.rp.seq)
    fp_rp_3p = three_prime_run(ps.fp.seq, ps.rp.seq)

    fp_gquad = 1.0 if has_gquad_seed_3p(ps.fp.seq) else 0.0
    rp_gquad = 1.0 if has_gquad_seed_3p(ps.rp.seq) else 0.0

    if fp_self > cfg.self_dimer_thresh_pct:
        ps.flags["fp_self_dimer"] = f"{fp_self:.1f}%"
    if rp_self > cfg.self_dimer_thresh_pct:
        ps.flags["rp_self_dimer"] = f"{rp_self:.1f}%"
    if fp_rp_cross > cfg.cross_dimer_thresh_pct:
        ps.flags["fp_rp_cross_dimer"] = f"{fp_rp_cross:.1f}%"
    if delta_tm > cfg.tm_delta_warn:
        ps.flags["delta_tm_high"] = f"{delta_tm:.1f}C"
    if fp_rp_3p >= 4:
        ps.flags["fp_rp_3p_run_ge4"] = str(fp_rp_3p)
    if fp_3p_self >= 4:
        ps.flags["fp_3p_self_ge4"] = str(fp_3p_self)
    if rp_3p_self >= 4:
        ps.flags["rp_3p_self_ge4"] = str(rp_3p_self)

    ps.metrics.update({
        "tm_fp_C": tm_fp,
        "tm_rp_C": tm_rp,
        "delta_tm_C": delta_tm,
        "fp_self_dimer_pct": fp_self,
        "rp_self_dimer_pct": rp_self,
        "fp_rp_cross_dimer_pct": fp_rp_cross,
        "fp_3p_self_run": float(fp_3p_self),
        "rp_3p_self_run": float(rp_3p_self),
        "fp_rp_3p_cross_run": float(fp_rp_3p),
        "fp_gquad_3p": fp_gquad,
        "rp_gquad_3p": rp_gquad,
        "fp_max_hpoly_run": float(max_homopolymer_run(ps.fp.seq)),
        "rp_max_hpoly_run": float(max_homopolymer_run(ps.rp.seq)),
        "fp_tail_gc": float(ps.fp.tail_gc_count),
        "rp_tail_gc": float(ps.rp.tail_gc_count),
    })

def annotate_context(ps: PairedSet, cfg: Config, target: str):
    amp_start = ps.fp.start
    amp_end = ps.rp.start + ps.rp.length  # exclusive
    amp_mid = amp_start + ps.amplicon_len / 2.0

    g_start = ps.guide.protospacer_start
    g_end = g_start + cfg.crrna_len
    g_mid = (g_start + g_end) / 2.0

    left_margin = g_start - amp_start
    right_margin = amp_end - g_end
    centered_abs_bp = abs(g_mid - amp_mid)

    if ps.guide.cas_type == "cas9":
        seed_seq = ps.guide.protospacer_seq[-8:]
    else:
        seed_seq = ps.guide.protospacer_seq[:8]
    seed_gc = gc_pct(seed_seq)
    seed_run = max_homopolymer_run(seed_seq)

    amp_seq = target[amp_start:amp_end]
    amp_gc = gc_pct(amp_seq)
    amp_run = max_homopolymer_run(amp_seq)

    pam_start = ps.guide.pam_start
    pam_len = len(ps.guide.pam_seq) if ps.guide.pam_seq else 0
    pam_end = pam_start + pam_len
    fp_span = (ps.fp.start, ps.fp.start + ps.fp.length)
    rp_span = (ps.rp.start, ps.rp.start + ps.rp.length)
    guide_span = (g_start, g_end)

    def spans_overlap(a: Tuple[int,int], b: Tuple[int,int]) -> bool:
        return not (a[1] <= b[0] or b[1] <= a[0])

    overlap_protospacer = (spans_overlap(fp_span, guide_span) or spans_overlap(rp_span, guide_span))
    overlap_pam = (pam_len > 0) and (spans_overlap(fp_span, (pam_start, pam_end)) or spans_overlap(rp_span, (pam_start, pam_end)))

    ps.metrics.update({
        "left_margin_bp": float(left_margin),
        "right_margin_bp": float(right_margin),
        "amplicon_center_offset_bp": float(centered_abs_bp),

        "seed_seq": seed_seq,
        "seed_len": 8.0,
        "seed_gc_pct": seed_gc,
        "seed_max_run": float(seed_run),
        "seed_has_homopolymer": 1.0 if seed_run >= 4 else 0.0,

        "amplicon_gc_pct": amp_gc,
        "amplicon_max_run": float(amp_run),

        "overlap_protospacer": 1.0 if overlap_protospacer else 0.0,
        "overlap_pam": 1.0 if overlap_pam else 0.0,
    })
    if overlap_protospacer:
        ps.flags["primer_overlaps_protospacer"] = "1"
    if overlap_pam:
        ps.flags["primer_overlaps_pam"] = "1"

def annotate_cross_with_guide(ps: PairedSet):
    gseq = ps.guide.protospacer_seq
    fp_g_pct = max_complementarity_percent(ps.fp.seq, gseq)
    rp_g_pct = max_complementarity_percent(ps.rp.seq, gseq)
    fp_g_run = three_prime_run(ps.fp.seq, gseq)
    rp_g_run = three_prime_run(ps.rp.seq, gseq)
    ps.metrics.update({
        "fp_guide_cross_pct": fp_g_pct,
        "rp_guide_cross_pct": rp_g_pct,
        "fp_guide_3p_run": float(fp_g_run),
        "rp_guide_3p_run": float(rp_g_run),
    })

def annotate_on_target_uniqueness(ps: PairedSet, target: str, cfg: Config):
    t = target.upper()
    k = min(cfg.on_target_k, ps.fp.length, ps.rp.length)

    rp_site_top = revcomp(ps.rp.seq)

    fp_hits = len(find_all_indices(t, ps.fp.seq)) + len(find_all_indices(t, revcomp(ps.fp.seq)))
    rp_hits = len(find_all_indices(t, rp_site_top)) + len(find_all_indices(t, revcomp(rp_site_top)))
    g_hits  = len(find_all_indices(t, ps.guide.protospacer_seq)) + len(find_all_indices(t, revcomp(ps.guide.protospacer_seq)))

    fp_addl_exact = max(0, fp_hits - 1)
    rp_addl_exact = max(0, rp_hits - 1)
    g_addl_exact  = max(0, g_hits - 1)

    fp_k = ps.fp.seq[-k:]
    rp_k = rp_site_top[-k:]
    fp_k_hits = len(find_all_indices(t, fp_k))
    rp_k_hits = len(find_all_indices(t, rp_k))

    fp_mm1 = count_approx_matches(t, ps.fp.seq, 1, extra_string=revcomp(ps.fp.seq))
    rp_mm1 = count_approx_matches(t, rp_site_top, 1, extra_string=revcomp(rp_site_top))
    g_mm1  = count_approx_matches(t, ps.guide.protospacer_seq, 1, extra_string=revcomp(ps.guide.protospacer_seq))

    fp_mm2 = count_approx_matches(t, ps.fp.seq, 2, extra_string=revcomp(ps.fp.seq))
    rp_mm2 = count_approx_matches(t, rp_site_top, 2, extra_string=revcomp(rp_site_top))
    g_mm2  = count_approx_matches(t, ps.guide.protospacer_seq, 2, extra_string=revcomp(ps.guide.protospacer_seq))

    fp_mm1_excl = max(0, fp_mm1 - 1)
    rp_mm1_excl = max(0, rp_mm1 - 1)
    g_mm1_excl  = max(0, g_mm1  - 1)
    fp_mm2_excl = max(0, fp_mm2 - 1)
    rp_mm2_excl = max(0, rp_mm2 - 1)
    g_mm2_excl  = max(0, g_mm2  - 1)

    fp_k_mm1 = count_approx_matches(t, fp_k, 1)
    rp_k_mm1 = count_approx_matches(t, rp_k, 1)

    min_amp, max_amp = cfg.amplicon_len_range
    rp_len = len(ps.rp.seq)
    fp_starts = find_all_indices(t, ps.fp.seq)
    rp_starts = find_all_indices(t, rp_site_top)
    amplicon_hits = 0
    for i in fp_starts:
        for j in rp_starts:
            if j <= i:
                continue
            amp_len = (j + rp_len) - i
            if min_amp <= amp_len <= max_amp:
                amplicon_hits += 1

    ps.metrics.update({
        "fp_full_hits_target_both": float(fp_hits),
        "rp_full_hits_target_both": float(rp_hits),
        "crrna_full_hits_target_both": float(g_hits),

        "fp_additional_exact_on_target": float(fp_addl_exact),
        "rp_additional_exact_on_target": float(rp_addl_exact),
        "crrna_additional_exact_on_target": float(g_addl_exact),

        "fp_full_hits_mm1_both": float(fp_mm1),
        "rp_full_hits_mm1_both": float(rp_mm1),
        "crrna_full_hits_mm1_both": float(g_mm1),

        "fp_full_hits_mm1_excl": float(fp_mm1_excl),
        "rp_full_hits_mm1_excl": float(rp_mm1_excl),
        "crrna_full_hits_mm1_excl": float(g_mm1_excl),

        "fp_full_hits_mm2_excl": float(fp_mm2_excl),
        "rp_full_hits_mm2_excl": float(rp_mm2_excl),
        "crrna_full_hits_mm2_excl": float(g_mm2_excl),

        "fp_3p_k": float(k),
        "fp_3p_k_hits_target": float(fp_k_hits),
        "rp_3p_k_hits_target": float(rp_k_hits),
        "fp_3p_k_hits_mm1": float(fp_k_mm1),
        "rp_3p_k_hits_mm1": float(rp_k_mm1),

        "amplicon_hits_on_target": float(amplicon_hits),
    })

    if fp_k_hits > 1:
        ps.flags["fp_3p_k_nonunique"] = f"{int(fp_k_hits)}"
    if rp_k_hits > 1:
        ps.flags["rp_3p_k_nonunique"] = f"{int(rp_k_hits)}"
    if amplicon_hits != 1:
        ps.flags["amplicon_not_unique"] = str(int(amplicon_hits))

# ---------- guide secondary structure ----------

def rnafold_cli_fold(seq_rna: str) -> Tuple[Optional[float], Optional[str]]:
    if which("RNAfold") is None:
        return (None, None)
    try:
        p = subprocess.run(
            ["RNAfold", "--noPS"],
            input=(seq_rna + "\n").encode("utf-8"),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
        out = p.stdout.decode("utf-8").strip().splitlines()
        if len(out) >= 2:
            struct_line = out[1]
            m = re.search(r"([().]+)\s+\(([-0-9.]+)\)", struct_line)
            if m:
                dot = m.group(1)
                mfe = float(m.group(2))
                return (mfe, dot)
    except Exception:
        pass
    return (None, None)

def viennarna_mfe_fold(seq_rna: str) -> Tuple[Optional[float], Optional[str]]:
    try:
        import RNA  # ViennaRNA Python module
    except Exception:
        return (None, None)
    try:
        fc = RNA.fold_compound(seq_rna)
        dot, mfe = fc.mfe()   # (dot-bracket, kcal/mol)
        return (float(mfe), str(dot))
    except Exception:
        return (None, None)

def rna_secondary_structure(seq_dna: str) -> Tuple[Optional[float], Optional[str]]:
    seq_rna = seq_dna.replace("T", "U")
    mfe, dot = viennarna_mfe_fold(seq_rna)
    if mfe is not None and dot is not None and len(dot) == len(seq_dna):
        return (mfe, dot)
    mfe, dot = rnafold_cli_fold(seq_rna)
    if mfe is not None and dot is not None and len(dot) == len(seq_dna):
        return (mfe, dot)
    return (None, None)

def annotate_guide_secondary_structure(ps: PairedSet, cfg: Config):
    g = ps.guide.protospacer_seq
    mfe, dot = rna_secondary_structure(g)
    seed_slice = slice(-8, None) if ps.guide.cas_type == "cas9" else slice(0, 8)

    if dot is not None and len(dot) == len(g):
        seed = dot[seed_slice]
        unpaired = seed.count(".")
        frac_unpaired = unpaired / max(1, len(seed))
        ps.metrics.update({
            "guide_mfe_kcal": float(mfe) if mfe is not None else 0.0,
            "guide_dotbracket": dot,
            "guide_seed_unpaired_frac": float(frac_unpaired),
        })
    else:
        seed_seq = (g[-8:] if ps.guide.cas_type == "cas9" else g[:8])
        best = 0
        rest = g[:-8] if ps.guide.cas_type != "cas9" else g[:-8]
        seed_rc = revcomp(seed_seq)
        for L in range(len(seed_seq), 0, -1):
            if seed_rc[:L] in rest:
                best = L
                break
        approx_unpaired_frac = max(0.0, 1.0 - best / max(1, len(seed_seq)))
        ps.metrics.update({
            "guide_mfe_kcal": 0.0,
            "guide_dotbracket": "",
            "guide_seed_unpaired_frac": float(approx_unpaired_frac),
        })

# =========================
# BLAST (optional)
# =========================

Hit = Dict[str, Union[float, int, str]]

class BlastCache:
    def __init__(self):
        self.cache: Dict[str, List[Hit]] = {}

    def get(self, seq: str):
        return self.cache.get(seq)

    def set(self, seq: str, val: List[Hit]):
        self.cache[seq] = val

def ensure_blast_db(background_fasta: str, workdir: str) -> str:
    if which("makeblastdb") is None:
        raise RuntimeError("makeblastdb not found in PATH.")
    db_prefix = os.path.join(workdir, "background_db")
    cmd = [
        "makeblastdb", "-in", background_fasta, "-dbtype", "nucl",
        "-parse_seqids", "-out", db_prefix
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return db_prefix

def blast_collect_hits(query_seq: str,
                       db_prefix: str,
                       cfg: Config,
                       workdir: str) -> List[Hit]:
    cachedir = os.path.join(workdir, "q")
    os.makedirs(cachedir, exist_ok=True)
    qfa = os.path.join(cachedir, "q.fa")
    with open(qfa, "w") as f:
        f.write(">q\n")
        f.write(query_seq + "\n")

    out_tsv = os.path.join(cachedir, "hits.tsv")

    cmd = [
        "blastn",
        "-query", qfa,
        "-db", db_prefix,
        "-word_size", str(cfg.blast_word_size),
        "-reward", "1", "-penalty", "-2",
        "-gapopen", "5", "-gapextend", "2",
        "-evalue", str(cfg.blast_evalue),
        "-strand", "both",
        "-task", cfg.blast_task,
        "-dust", "no",
        "-soft_masking", "false",
        "-outfmt", "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore sstrand",
        "-out", out_tsv
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    hits: List[Hit] = []
    qlen = int(len(query_seq))
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 12:
                    continue
                pident   = float(parts[2])
                aln_len  = int(parts[3])
                mismatch = int(parts[4])
                qstart   = int(parts[5])
                qend     = int(parts[6])
                sstart   = int(parts[7])
                send     = int(parts[8])
                sseqid   = parts[1]
                sstrand  = parts[11]
                score    = (pident/100.0) * (aln_len / max(1.0, float(qlen)))
                is_exact_full = (mismatch == 0 and aln_len == qlen and ((qstart == 1 and qend == qlen) or (qstart == qlen and qend == 1)))
                hits.append({
                    "sseqid": sseqid, "pident": pident, "aln_len": aln_len, "mismatch": mismatch,
                    "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
                    "sstrand": sstrand, "score": score, "is_exact_full": 1 if is_exact_full else 0
                })
    return hits

def _top_any_and_imperfect(hits: List[Hit], qlen: int):
    if not hits:
        return None, None, 0, ""
    hits_sorted = sorted(hits, key=lambda h: h["score"], reverse=True)
    top_any = hits_sorted[0]
    n_exact = sum(1 for h in hits_sorted if int(h["is_exact_full"]) == 1)
    example_exact_id = ""
    for h in hits_sorted:
        if int(h["is_exact_full"]) == 1:
            example_exact_id = str(h["sseqid"])
            break
    imperfect = None
    for h in hits_sorted:
        if int(h["mismatch"]) > 0 or int(h["aln_len"]) != qlen:
            imperfect = h
            break
    return top_any, imperfect, n_exact, example_exact_id

def annotate_specificity_with_blast(pairs: List[PairedSet], cfg: Config, background_fasta: str):
    if which("blastn") is None:
        raise RuntimeError("blastn not found in PATH.")
    with tempfile.TemporaryDirectory(prefix="primedrpa_blast_") as td:
        db_prefix = ensure_blast_db(background_fasta, td)
        cache = BlastCache()

        def get_hits(seq: str):
            cached = cache.get(seq)
            if cached is not None:
                return cached
            h = blast_collect_hits(seq, db_prefix, cfg, td)
            cache.set(seq, h)
            return h

        for ps in pairs:
            fp_hits = get_hits(ps.fp.seq)
            rp_hits = get_hits(ps.rp.seq)
            g_hits  = get_hits(ps.guide.protospacer_seq)

            qlen_fp = len(ps.fp.seq)
            qlen_rp = len(ps.rp.seq)
            qlen_g  = len(ps.guide.protospacer_seq)

            fp_any, fp_imp, fp_exact_n, fp_exact_id = _top_any_and_imperfect(fp_hits, qlen_fp)
            rp_any, rp_imp, rp_exact_n, rp_exact_id = _top_any_and_imperfect(rp_hits, qlen_rp)
            g_any,  g_imp,  g_exact_n,  g_exact_id  = _top_any_and_imperfect(g_hits,  qlen_g)

            def _save(prefix: str, top, qlen: int, exact_n: int, exact_id: str):
                if top is None:
                    ps.metrics[f"{prefix}_offtarget_score"] = 0.0
                    ps.metrics[f"{prefix}_offtarget_hit"] = ""
                    ps.metrics[f"{prefix}_offtarget_mm"] = qlen
                    ps.metrics[f"{prefix}_offtarget_aln"] = 0
                else:
                    ps.metrics[f"{prefix}_offtarget_score"] = float(top["score"])
                    ps.metrics[f"{prefix}_offtarget_hit"] = str(top["sseqid"])
                    ps.metrics[f"{prefix}_offtarget_mm"] = int(top["mismatch"])
                    ps.metrics[f"{prefix}_offtarget_aln"] = int(top["aln_len"])
                ps.metrics[f"{prefix}_has_exact_full_hit"] = 1 if exact_n > 0 else 0
                ps.metrics[f"{prefix}_exact_full_hit_count"] = int(exact_n)
                ps.metrics[f"{prefix}_exact_full_example_sseqid"] = exact_id

            def _save_imp(prefix: str, top_imp, top_any, qlen: int):
                src = top_imp if top_imp is not None else top_any
                tag = "imperfect" if top_imp is not None else "any"
                if src is None:
                    ps.metrics[f"{prefix}_offtarget_score_imp"] = 0.0
                    ps.metrics[f"{prefix}_offtarget_mm_imp"] = qlen
                    ps.metrics[f"{prefix}_offtarget_from"] = "none"
                else:
                    ps.metrics[f"{prefix}_offtarget_score_imp"] = float(src["score"])
                    ps.metrics[f"{prefix}_offtarget_mm_imp"] = int(src["mismatch"])
                    ps.metrics[f"{prefix}_offtarget_from"] = tag

            _save("fp", fp_any, qlen_fp, fp_exact_n, fp_exact_id); _save_imp("fp", fp_imp, fp_any, qlen_fp)
            _save("rp", rp_any, qlen_rp, rp_exact_n, rp_exact_id); _save_imp("rp", rp_imp, rp_any, qlen_rp)
            _save("crrna", g_any, qlen_g, g_exact_n, g_exact_id); _save_imp("crrna", g_imp, g_any, qlen_g)

# ---------- conservation annotator ----------

def _ungapped_to_aln_cols(ref_map: List[int], start_ungap: int, end_ungap_excl: int) -> List[int]:
    cols = []
    for i in range(start_ungap, end_ungap_excl):
        if 0 <= i < len(ref_map):
            cols.append(ref_map[i])
    return cols

def annotate_conservation(ps: PairedSet,
                          msa_info: Optional[Dict[str, List[float]]],
                          cfg: Config):
    if not msa_info:
        ps.metrics.update({
            "guide_conservation": 0.5,
            "fp_conservation": 0.5,
            "rp_conservation": 0.5,
            "amplicon_conservation": 0.5,
        })
        return
    cons_cols = msa_info["cons_cols"]
    ref_map = msa_info["ref_map"]

    g_start = ps.guide.protospacer_start
    g_end = g_start + cfg.crrna_len
    fp_start = ps.fp.start; fp_end = fp_start + ps.fp.length
    rp_start = ps.rp.start; rp_end = rp_start + ps.rp.length
    amp_start = ps.fp.start; amp_end = ps.rp.start + ps.rp.length

    def region_cons(start_u, end_u):
        cols = _ungapped_to_aln_cols(ref_map, start_u, end_u)
        vals = [cons_cols[c] for c in cols if 0 <= c < len(cons_cols)]
        return _mean(vals) if vals else 0.5

    ps.metrics.update({
        "guide_conservation": float(region_cons(g_start, g_end)),
        "fp_conservation": float(region_cons(fp_start, fp_end)),
        "rp_conservation": float(region_cons(rp_start, rp_end)),
        "amplicon_conservation": float(region_cons(amp_start, amp_end)),
    })

# ---------- PAM subscore ----------

def pam_subscore(ps: PairedSet) -> float:
    overlap = int(ps.metrics.get("overlap_pam", 0))
    if overlap:
        return 0.0
    if ps.guide.cas_type == "cas9":
        return 1.0 if len(ps.guide.pam_seq) == 3 and ps.guide.pam_seq[1:] == "GG" else 0.7
    if ps.guide.cas_type == "cas12a":
        return 1.0 if len(ps.guide.pam_seq) == 4 and ps.guide.pam_seq[:3] == "TTT" else 0.7
    return 1.0

# ---------- specificity fallback (no background) ----------

def target_specificity_fallback(ps: PairedSet) -> float:
    m = ps.metrics
    amp_ok = _less_is_better(abs(m.get("amplicon_hits_on_target", 1.0) - 1.0), good=0.0, bad=2.0)
    extra_exact = min(
        _less_is_better(m.get("fp_additional_exact_on_target", 0.0), 0.0, 1.0),
        _less_is_better(m.get("rp_additional_exact_on_target", 0.0), 0.0, 1.0),
        _less_is_better(m.get("crrna_additional_exact_on_target", 0.0), 0.0, 1.0),
    )
    approx1 = min(
        _less_is_better(m.get("fp_full_hits_mm1_excl", 0.0), 0.0, 3.0),
        _less_is_better(m.get("rp_full_hits_mm1_excl", 0.0), 0.0, 3.0),
        _less_is_better(m.get("crrna_full_hits_mm1_excl", 0.0), 0.0, 3.0),
    )
    approx2 = min(
        _less_is_better(m.get("fp_full_hits_mm2_excl", 0.0), 0.0, 6.0),
        _less_is_better(m.get("rp_full_hits_mm2_excl", 0.0), 0.0, 6.0),
        _less_is_better(m.get("crrna_full_hits_mm2_excl", 0.0), 0.0, 6.0),
    )
    seed_nonuniq = min(
        _less_is_better(max(0.0, m.get("fp_3p_k_hits_target", 1.0) - 1.0), 0.0, 2.0),
        _less_is_better(max(0.0, m.get("rp_3p_k_hits_target", 1.0) - 1.0), 0.0, 2.0),
    )
    return min(amp_ok, extra_exact, approx1, approx2, seed_nonuniq)

# ---------- composite scoring ----------

def compute_scores(ps: PairedSet, cfg: Config, background_used: bool, weights: Optional[Dict[str, float]] = None):
    m = ps.metrics
    g = ps.guide

    # Geometry
    geo = min(
        _less_is_better(abs(m.get("amplicon_center_offset_bp", 0.0)), good=0.0, bad=30.0),
        _more_is_better(m.get("left_margin_bp", 0.0), good=cfg.guide_margin_bp, cap=cfg.guide_margin_bp*2),
        _more_is_better(m.get("right_margin_bp", 0.0), good=cfg.guide_margin_bp, cap=cfg.guide_margin_bp*2),
    )

    # Primers thermo
    tm_fp = m.get("tm_fp_C", 0.0)
    tm_rp = m.get("tm_rp_C", 0.0)
    delta_tm = abs(m.get("delta_tm_C", 0.0))
    tm_target = 63.0
    thermo = min(
        _tri_range(tm_fp, ideal=tm_target, tol=6.0, hard_lo=50.0, hard_hi=75.0),
        _tri_range(tm_rp, ideal=tm_target, tol=6.0, hard_lo=50.0, hard_hi=75.0),
        _less_is_better(delta_tm, good=0.0, bad=5.0)
    )

    # Structure
    struct_ok = min(
        _less_is_better(m.get("fp_self_dimer_pct", 0.0), good=0.0, bad=cfg.self_dimer_thresh_pct),
        _less_is_better(m.get("rp_self_dimer_pct", 0.0), good=0.0, bad=cfg.self_dimer_thresh_pct),
        _less_is_better(m.get("fp_rp_cross_dimer_pct", 0.0), good=0.0, bad=cfg.cross_dimer_thresh_pct),
        _less_is_better(m.get("fp_3p_self_run", 0.0), good=0.0, bad=4.0),
        _less_is_better(m.get("rp_3p_self_run", 0.0), good=0.0, bad=4.0),
        _less_is_better(m.get("fp_rp_3p_cross_run", 0.0), good=0.0, bad=3.0),
        1.0 - min(1.0, (m.get("fp_gquad_3p", 0.0) + m.get("rp_gquad_3p", 0.0)) * 0.75)
    )

    # Context
    seed_gc = m.get("seed_gc_pct", 50.0)
    seed_run = m.get("seed_max_run", 0.0)
    amp_gc = m.get("amplicon_gc_pct", 50.0)
    context = min(
        _band(seed_gc, 40.0, 60.0, soft_margin=10.0),
        _band(amp_gc, 35.0, 65.0, soft_margin=10.0),
        _less_is_better(seed_run, good=0.0, bad=4.0)
    )

    # Uniqueness
    uniq = min(
        _less_is_better(m.get("fp_3p_k_hits_target", 0.0), good=1.0, bad=3.0),
        _less_is_better(m.get("rp_3p_k_hits_target", 0.0), good=1.0, bad=3.0),
        _less_is_better(m.get("fp_additional_exact_on_target", 0.0), good=0.0, bad=2.0),
        _less_is_better(m.get("rp_additional_exact_on_target", 0.0), good=0.0, bad=2.0),
        _less_is_better(m.get("crrna_additional_exact_on_target", 0.0), good=0.0, bad=2.0),
        _less_is_better(m.get("fp_full_hits_mm1_excl", 0.0), good=0.0, bad=4.0),
        _less_is_better(m.get("rp_full_hits_mm1_excl", 0.0), good=0.0, bad=4.0),
        _less_is_better(m.get("crrna_full_hits_mm1_excl", 0.0), good=0.0, bad=4.0),
        _less_is_better(m.get("fp_full_hits_mm2_excl", 0.0), good=0.0, bad=8.0),
        _less_is_better(m.get("rp_full_hits_mm2_excl", 0.0), good=0.0, bad=8.0),
        _less_is_better(m.get("crrna_full_hits_mm2_excl", 0.0), good=0.0, bad=8.0),
        _less_is_better(abs(m.get("amplicon_hits_on_target", 1.0) - 1.0), good=0.0, bad=2.0)
    )

    # Specificity: background BLAST or target fallback
    if background_used:
        spec = min(
            _less_is_better(m.get("fp_offtarget_score_imp", 0.0), good=0.0, bad=cfg.offtarget_max_allow),
            _less_is_better(m.get("rp_offtarget_score_imp", 0.0), good=0.0, bad=cfg.offtarget_max_allow),
            _less_is_better(m.get("crrna_offtarget_score_imp", 0.0), good=0.0, bad=cfg.offtarget_max_allow),
            _more_is_better(m.get("fp_offtarget_mm_imp", ps.fp.length), good=ps.fp.length, cap=ps.fp.length),
            _more_is_better(m.get("rp_offtarget_mm_imp", ps.rp.length), good=ps.rp.length, cap=ps.rp.length),
            _more_is_better(m.get("crrna_offtarget_mm_imp", len(g.protospacer_seq)), good=len(g.protospacer_seq), cap=len(g.protospacer_seq)),
        )
        for pfx in ("fp", "rp", "crrna"):
            if m.get(f"{pfx}_has_exact_full_hit", 0):
                spec = min(spec, 0.0)
    else:
        spec = target_specificity_fallback(ps)

    # Guide quality (GC + seed accessibility)
    guide_gc = g.gc_pct
    guide_struct = min(
        _band(guide_gc, 35.0, 60.0, soft_margin=10.0),
        _more_is_better(m.get("guide_seed_unpaired_frac", 0.0), good=0.7, cap=1.0)
    )

    # Conservation
    cons = _mean([
        m.get("guide_conservation", 0.5),
        m.get("fp_conservation", 0.5),
        m.get("rp_conservation", 0.5),
        m.get("amplicon_conservation", 0.5),
    ])

    # PAM/site context
    pam = pam_subscore(ps)

    # Default weights
    W = {"uniq": 0.20,"thermo": 0.16,"spec": 0.14,"struct": 0.12,"geo": 0.10,"guide": 0.10,"ctx": 0.07,"cons": 0.06,"pam": 0.05}
    if weights:
        for k in list(W.keys()):
            if k in weights:
                W[k] = float(weights[k])
        total = sum(W.values()) or 1.0
        for k in W:
            W[k] = W[k]/total

    score = (
        W["uniq"]*uniq + W["thermo"]*thermo + W["spec"]*spec + W["struct"]*struct_ok +
        W["geo"]*geo + W["guide"]*guide_struct + W["ctx"]*context + W["cons"]*cons + W["pam"]*pam
    )

    ps.metrics.update({
        "score_overall": float(score),
        "score_uniqueness": float(uniq),
        "score_primers_thermo": float(thermo),
        "score_specificity": float(spec),
        "score_structure_penalty": float(struct_ok),
        "score_geometry": float(geo),
        "score_guide": float(guide_struct),
        "score_context": float(context),
        "score_conservation": float(cons),
        "score_pam": float(pam),
    })

# =========================
# CSV writer (sorted; includes rank)
# =========================

def write_pairs_csv(pairs: List[PairedSet], out_path: str, target_name: str, background_used: bool):
    fields = [
        "rank",
        # identifiers
        "target", "cas_type", "guide_idx", "guide_strand",
        "guide_pam", "guide_pam_start", "guide_start", "guide_len",
        "guide_gc_pct", "guide_seq",
        # primers & geometry
        "fp_start", "fp_len", "fp_gc_pct", "fp_seq", "fp_tail_gc", "fp_max_hpoly_run",
        "rp_start", "rp_len", "rp_gc_pct", "rp_seq", "rp_tail_gc", "rp_max_hpoly_run",
        "amplicon_len",
        "left_margin_bp", "right_margin_bp", "amplicon_center_offset_bp",
        # thermo/structure
        "tm_fp_C", "tm_rp_C", "delta_tm_C",
        "fp_self_dimer_pct", "rp_self_dimer_pct", "fp_rp_cross_dimer_pct",
        "fp_3p_self_run", "rp_3p_self_run", "fp_rp_3p_cross_run",
        "fp_gquad_3p", "rp_gquad_3p",
        # seed/context
        "seed_seq", "seed_len", "seed_gc_pct", "seed_max_run", "seed_has_homopolymer",
        "amplicon_gc_pct", "amplicon_max_run",
        "overlap_protospacer", "overlap_pam",
        # cross-talk with guide
        "fp_guide_cross_pct", "rp_guide_cross_pct", "fp_guide_3p_run", "rp_guide_3p_run",
        # on-target uniqueness & approx (target-only)
        "fp_full_hits_target_both", "rp_full_hits_target_both", "crrna_full_hits_target_both",
        "fp_additional_exact_on_target", "rp_additional_exact_on_target", "crrna_additional_exact_on_target",
        "fp_full_hits_mm1_both", "rp_full_hits_mm1_both", "crrna_full_hits_mm1_both",
        "fp_full_hits_mm1_excl", "rp_full_hits_mm1_excl", "crrna_full_hits_mm1_excl",
        "fp_full_hits_mm2_excl", "rp_full_hits_mm2_excl", "crrna_full_hits_mm2_excl",
        "fp_3p_k", "fp_3p_k_hits_target", "rp_3p_k_hits_target",
        "fp_3p_k_hits_mm1", "rp_3p_k_hits_mm1",
        "amplicon_hits_on_target",
        # background usage marker
        "background_used",
        # mismatch-aware individual off-targets vs BACKGROUND
        "fp_offtarget_score", "fp_offtarget_hit", "fp_offtarget_mm", "fp_offtarget_aln",
        "fp_offtarget_score_imp", "fp_offtarget_mm_imp", "fp_offtarget_from",
        "fp_has_exact_full_hit", "fp_exact_full_hit_count", "fp_exact_full_example_sseqid",
        "rp_offtarget_score", "rp_offtarget_hit", "rp_offtarget_mm", "rp_offtarget_aln",
        "rp_offtarget_score_imp", "rp_offtarget_mm_imp", "rp_offtarget_from",
        "rp_has_exact_full_hit", "rp_exact_full_hit_count", "rp_exact_full_example_sseqid",
        "crrna_offtarget_score", "crrna_offtarget_hit", "crrna_offtarget_mm", "crrna_offtarget_aln",
        "crrna_offtarget_score_imp", "crrna_offtarget_mm_imp", "crrna_offtarget_from",
        "crrna_has_exact_full_hit", "crrna_exact_full_hit_count", "crrna_exact_full_example_sseqid",
        # guide secondary structure
        "guide_mfe_kcal", "guide_seed_unpaired_frac", "guide_dotbracket",
        # conservation
        "guide_conservation", "fp_conservation", "rp_conservation", "amplicon_conservation",
        # scores (objective, 0..1)
        "score_overall", "score_uniqueness", "score_primers_thermo", "score_specificity",
        "score_structure_penalty", "score_geometry", "score_guide", "score_context", "score_conservation", "score_pam",
        # flags
        "flags"
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for ps in pairs:
            g = ps.guide
            row = {
                "rank": int(ps.metrics.get("rank", 0)),

                "target": target_name,
                "cas_type": g.cas_type,
                "guide_idx": ps.guide_index,
                "guide_strand": g.strand,
                "guide_pam": g.pam_seq,
                "guide_pam_start": g.pam_start,
                "guide_start": g.protospacer_start,
                "guide_len": len(g.protospacer_seq),
                "guide_gc_pct": f"{g.gc_pct:.1f}",
                "guide_seq": g.protospacer_seq,

                "fp_start": ps.fp.start,
                "fp_len": ps.fp.length,
                "fp_gc_pct": f"{ps.fp.gc_pct:.1f}",
                "fp_seq": ps.fp.seq,
                "fp_tail_gc": int(ps.metrics.get("fp_tail_gc", 0)),
                "fp_max_hpoly_run": int(ps.metrics.get("fp_max_hpoly_run", 0)),

                "rp_start": ps.rp.start,
                "rp_len": ps.rp.length,
                "rp_gc_pct": f"{ps.rp.gc_pct:.1f}",
                "rp_seq": ps.rp.seq,
                "rp_tail_gc": int(ps.metrics.get("rp_tail_gc", 0)),
                "rp_max_hpoly_run": int(ps.metrics.get("rp_max_hpoly_run", 0)),

                "amplicon_len": ps.amplicon_len,
                "left_margin_bp": f"{ps.metrics.get('left_margin_bp', 0):.1f}",
                "right_margin_bp": f"{ps.metrics.get('right_margin_bp', 0):.1f}",
                "amplicon_center_offset_bp": f"{ps.metrics.get('amplicon_center_offset_bp', 0):.1f}",

                "tm_fp_C": f"{ps.metrics.get('tm_fp_C', 0):.1f}",
                "tm_rp_C": f"{ps.metrics.get('tm_rp_C', 0):.1f}",
                "delta_tm_C": f"{ps.metrics.get('delta_tm_C', 0):.1f}",
                "fp_self_dimer_pct": f"{ps.metrics.get('fp_self_dimer_pct', 0):.1f}",
                "rp_self_dimer_pct": f"{ps.metrics.get('rp_self_dimer_pct', 0):.1f}",
                "fp_rp_cross_dimer_pct": f"{ps.metrics.get('fp_rp_cross_dimer_pct', 0):.1f}",
                "fp_3p_self_run": int(ps.metrics.get('fp_3p_self_run', 0)),
                "rp_3p_self_run": int(ps.metrics.get('rp_3p_self_run', 0)),
                "fp_rp_3p_cross_run": int(ps.metrics.get('fp_rp_3p_cross_run', 0)),
                "fp_gquad_3p": int(ps.metrics.get('fp_gquad_3p', 0)),
                "rp_gquad_3p": int(ps.metrics.get('rp_gquad_3p', 0)),

                "seed_seq": ps.metrics.get("seed_seq", ""),
                "seed_len": int(ps.metrics.get("seed_len", 0)),
                "seed_gc_pct": f"{ps.metrics.get('seed_gc_pct', 0):.1f}",
                "seed_max_run": int(ps.metrics.get("seed_max_run", 0)),
                "seed_has_homopolymer": int(ps.metrics.get('seed_has_homopolymer', 0)),

                "amplicon_gc_pct": f"{ps.metrics.get('amplicon_gc_pct', 0):.1f}",
                "amplicon_max_run": int(ps.metrics.get("amplicon_max_run", 0)),

                "overlap_protospacer": int(ps.metrics.get("overlap_protospacer", 0)),
                "overlap_pam": int(ps.metrics.get("overlap_pam", 0)),

                "fp_guide_cross_pct": f"{ps.metrics.get('fp_guide_cross_pct', 0):.1f}",
                "rp_guide_cross_pct": f"{ps.metrics.get('rp_guide_cross_pct', 0):.1f}",
                "fp_guide_3p_run": int(ps.metrics.get('fp_guide_3p_run', 0)),
                "rp_guide_3p_run": int(ps.metrics.get('rp_guide_3p_run', 0)),

                "fp_full_hits_target_both": int(ps.metrics.get('fp_full_hits_target_both', 0)),
                "rp_full_hits_target_both": int(ps.metrics.get('rp_full_hits_target_both', 0)),
                "crrna_full_hits_target_both": int(ps.metrics.get('crrna_full_hits_target_both', 0)),

                "fp_additional_exact_on_target": int(ps.metrics.get('fp_additional_exact_on_target', 0)),
                "rp_additional_exact_on_target": int(ps.metrics.get('rp_additional_exact_on_target', 0)),
                "crrna_additional_exact_on_target": int(ps.metrics.get('crrna_additional_exact_on_target', 0)),

                "fp_full_hits_mm1_both": int(ps.metrics.get('fp_full_hits_mm1_both', 0)),
                "rp_full_hits_mm1_both": int(ps.metrics.get('rp_full_hits_mm1_both', 0)),
                "crrna_full_hits_mm1_both": int(ps.metrics.get('crrna_full_hits_mm1_both', 0)),

                "fp_full_hits_mm1_excl": int(ps.metrics.get('fp_full_hits_mm1_excl', 0)),
                "rp_full_hits_mm1_excl": int(ps.metrics.get('rp_full_hits_mm1_excl', 0)),
                "crrna_full_hits_mm1_excl": int(ps.metrics.get('crrna_full_hits_mm1_excl', 0)),

                "fp_full_hits_mm2_excl": int(ps.metrics.get('fp_full_hits_mm2_excl', 0)),
                "rp_full_hits_mm2_excl": int(ps.metrics.get('rp_full_hits_mm2_excl', 0)),
                "crrna_full_hits_mm2_excl": int(ps.metrics.get('crrna_full_hits_mm2_excl', 0)),

                "fp_3p_k": int(ps.metrics.get("fp_3p_k", 0)),
                "fp_3p_k_hits_target": int(ps.metrics.get('fp_3p_k_hits_target', 0)),
                "rp_3p_k_hits_target": int(ps.metrics.get('rp_3p_k_hits_target', 0)),
                "fp_3p_k_hits_mm1": int(ps.metrics.get('fp_3p_k_hits_mm1', 0)),
                "rp_3p_k_hits_mm1": int(ps.metrics.get('rp_3p_k_hits_mm1', 0)),

                "amplicon_hits_on_target": int(ps.metrics.get('amplicon_hits_on_target', 0)),

                "background_used": 1 if background_used else 0,

                "fp_offtarget_score": f"{ps.metrics.get('fp_offtarget_score', 0.0):.3f}",
                "fp_offtarget_hit": ps.metrics.get("fp_offtarget_hit",""),
                "fp_offtarget_mm": int(ps.metrics.get("fp_offtarget_mm", 0)),
                "fp_offtarget_aln": int(ps.metrics.get("fp_offtarget_aln", 0)),
                "fp_offtarget_score_imp": f"{ps.metrics.get('fp_offtarget_score_imp', 0.0):.3f}",
                "fp_offtarget_mm_imp": int(ps.metrics.get("fp_offtarget_mm_imp", 0)),
                "fp_offtarget_from": ps.metrics.get("fp_offtarget_from",""),
                "fp_has_exact_full_hit": int(ps.metrics.get("fp_has_exact_full_hit", 0)),
                "fp_exact_full_hit_count": int(ps.metrics.get("fp_exact_full_hit_count", 0)),
                "fp_exact_full_example_sseqid": ps.metrics.get("fp_exact_full_example_sseqid",""),

                "rp_offtarget_score": f"{ps.metrics.get('rp_offtarget_score', 0.0):.3f}",
                "rp_offtarget_hit": ps.metrics.get("rp_offtarget_hit",""),
                "rp_offtarget_mm": int(ps.metrics.get("rp_offtarget_mm", 0)),
                "rp_offtarget_aln": int(ps.metrics.get("rp_offtarget_aln", 0)),
                "rp_offtarget_score_imp": f"{ps.metrics.get('rp_offtarget_score_imp', 0.0):.3f}",
                "rp_offtarget_mm_imp": int(ps.metrics.get("rp_offtarget_mm_imp", 0)),
                "rp_offtarget_from": ps.metrics.get("rp_offtarget_from",""),
                "rp_has_exact_full_hit": int(ps.metrics.get("rp_has_exact_full_hit", 0)),
                "rp_exact_full_hit_count": int(ps.metrics.get("rp_exact_full_hit_count", 0)),
                "rp_exact_full_example_sseqid": ps.metrics.get("rp_exact_full_example_sseqid",""),

                "crrna_offtarget_score": f"{ps.metrics.get('crrna_offtarget_score', 0.0):.3f}",
                "crrna_offtarget_hit": ps.metrics.get("crrna_offtarget_hit",""),
                "crrna_offtarget_mm": int(ps.metrics.get("crrna_offtarget_mm", 0)),
                "crrna_offtarget_aln": int(ps.metrics.get("crrna_offtarget_aln", 0)),
                "crrna_offtarget_score_imp": f"{ps.metrics.get('crrna_offtarget_score_imp', 0.0):.3f}",
                "crrna_offtarget_mm_imp": int(ps.metrics.get("crrna_offtarget_mm_imp", 0)),
                "crrna_offtarget_from": ps.metrics.get("crrna_offtarget_from",""),
                "crrna_has_exact_full_hit": int(ps.metrics.get("crrna_has_exact_full_hit", 0)),
                "crrna_exact_full_hit_count": int(ps.metrics.get("crrna_exact_full_hit_count", 0)),
                "crrna_exact_full_example_sseqid": ps.metrics.get("crrna_exact_full_example_sseqid",""),

                "guide_mfe_kcal": f"{ps.metrics.get('guide_mfe_kcal', 0):.2f}",
                "guide_seed_unpaired_frac": f"{ps.metrics.get('guide_seed_unpaired_frac', 0):.3f}",
                "guide_dotbracket": ps.metrics.get("guide_dotbracket", ""),

                "guide_conservation": f"{ps.metrics.get('guide_conservation', 0.5):.3f}",
                "fp_conservation": f"{ps.metrics.get('fp_conservation', 0.5):.3f}",
                "rp_conservation": f"{ps.metrics.get('rp_conservation', 0.5):.3f}",
                "amplicon_conservation": f"{ps.metrics.get('amplicon_conservation', 0.5):.3f}",

                "score_overall": f"{ps.metrics.get('score_overall', 0.0):.4f}",
                "score_uniqueness": f"{ps.metrics.get('score_uniqueness', 0.0):.4f}",
                "score_primers_thermo": f"{ps.metrics.get('score_primers_thermo', 0.0):.4f}",
                "score_specificity": f"{ps.metrics.get('score_specificity', 0.0):.4f}",
                "score_structure_penalty": f"{ps.metrics.get('score_structure_penalty', 0.0):.4f}",
                "score_geometry": f"{ps.metrics.get('score_geometry', 0.0):.4f}",
                "score_guide": f"{ps.metrics.get('score_guide', 0.0):.4f}",
                "score_context": f"{ps.metrics.get('score_context', 0.0):.4f}",
                "score_conservation": f"{ps.metrics.get('score_conservation', 0.0):.4f}",
                "score_pam": f"{ps.metrics.get('score_pam', 0.0):.4f}",

                "flags": ";".join([f"{k}={v}" for k, v in ps.flags.items()]) if ps.flags else ""
            }
            w.writerow(row)
    print(f"[WRITE] CSV → {out_path}  ({len(pairs)} rows, sorted by score)")

# =========================
# Main
# =========================

if __name__ == "__main__":
    try:
        cfg, tname, tseq, bg, out_csv, msa_info, weights = step1_get_config_and_sequences()

        # Step 2
        guides = step2_discover_crrnas(tseq, cfg)
        if not guides:
            print("\n[NOTE] No crRNA candidates passed light filtering.")
            sys.exit(0)

        # Step 3
        fps = candidate_primers_plus_strand(tseq, cfg)
        rps = candidate_primers_minus_strand(tseq, cfg)
        print(
            f"\n[PRIMERS] {len(fps)} forward; {len(rps)} reverse candidates passed cheap filters "
            f"(GC {cfg.gc_primer_range[0]}–{cfg.gc_primer_range[1]}%, homopolymer<{cfg.max_homopolymer}, 3' clamp=1 GC)."
        )
        if not fps or not rps:
            print("[NOTE] No primer candidates after cheap filters.")
            sys.exit(0)

        # Step 4
        pairs = pair_primers_around_guides(tseq, guides, fps, rps, cfg)
        print(
            f"[PAIRING] {len(pairs)} FP/RP pairs span guides with margin {cfg.guide_margin_bp} bp and amplicon "
            f"{cfg.amplicon_len_range[0]}–{cfg.amplicon_len_range[1]} bp."
        )
        if not pairs:
            print("[NOTE] No pairs satisfied geometry/amplicon constraints.")
            sys.exit(0)

        # Steps 5–6: features
        for ps in pairs:
            annotate_structure_and_tm(ps, cfg)
            annotate_context(ps, cfg, tseq)
            annotate_cross_with_guide(ps)
            annotate_on_target_uniqueness(ps, tseq, cfg)
            annotate_guide_secondary_structure(ps, cfg)
            annotate_conservation(ps, msa_info, cfg)

        # Step 7: Optional BLAST specificity (per-component only)
        background_used = bool(bg)
        if background_used:
            try:
                print("[BLAST] Building DB and checking off-targets for FP, RP, and crRNA (per-component, mismatch-aware)...")
                annotate_specificity_with_blast(pairs, cfg, bg)
            except Exception as e:
                print(f"[WARN] Skipping BLAST specificity due to error: {e}")
                background_used = False  # fall back to target-based specificity scoring

        # Scoring (objective, fully weighted)
        for ps in pairs:
            compute_scores(ps, cfg, background_used=background_used, weights=weights)

        # ---------- NEW: RANKING ----------
        # Sort by composite score (desc), assign ranks, and show a top-K preview.
        pairs_sorted = sorted(pairs, key=lambda p: p.metrics.get("score_overall", 0.0), reverse=True)
        for rank_idx, ps in enumerate(pairs_sorted, start=1):
            ps.metrics["rank"] = rank_idx

        show_n = min(cfg.top_k, len(pairs_sorted))
        print(f"\n[RANKED] Showing top {show_n} of {len(pairs_sorted)} (highest → lowest score):")
        for ps in pairs_sorted[:show_n]:
            g = ps.guide
            print(
                f"\n#{int(ps.metrics['rank'])}  score={ps.metrics['score_overall']:.3f}   "
                f"[uniq {ps.metrics['score_uniqueness']:.2f} | thermo {ps.metrics['score_primers_thermo']:.2f} | "
                f"spec {ps.metrics['score_specificity']:.2f} | struct {ps.metrics['score_structure_penalty']:.2f} | "
                f"geo {ps.metrics['score_geometry']:.2f} | guide {ps.metrics['score_guide']:.2f} | "
                f"ctx {ps.metrics['score_context']:.2f} | cons {ps.metrics['score_conservation']:.2f} | pam {ps.metrics['score_pam']:.2f}]"
            )
            pam_info = f"PAM={g.pam_seq} at {g.pam_start}" if g.pam_seq else "PAM=NA"
            print(
                f"  GuideIdx={ps.guide_index} {g.cas_type}{g.strand}  {pam_info}  prot@{g.protospacer_start} "
                f"len={len(g.protospacer_seq)} GC={g.gc_pct:.1f}%  seed={ps.metrics['seed_seq']}"
            )
            print(
                f"  FP  start={ps.fp.start} len={ps.fp.length} GC={ps.fp.gc_pct:.1f}%  "
                f"Tm={ps.metrics['tm_fp_C']:.1f}°C  self%={ps.metrics['fp_self_dimer_pct']:.1f}  seq={ps.fp.seq}"
            )
            print(
                f"  RP  start={ps.rp.start} len={ps.rp.length} GC={ps.rp.gc_pct:.1f}%  "
                f"Tm={ps.metrics['tm_rp_C']:.1f}°C  self%={ps.metrics['rp_self_dimer_pct']:.1f}  seq={ps.rp.seq}"
            )
            print(
                f"  Cross-dimer%={ps.metrics['fp_rp_cross_dimer_pct']:.1f}   "
                f"ΔTm={ps.metrics['delta_tm_C']:.1f}°C   3'cross_run={int(ps.metrics['fp_rp_3p_cross_run'])}"
            )
            lm = ps.metrics['left_margin_bp']; rm = ps.metrics['right_margin_bp']
            print(
                f"  Amplicon={ps.amplicon_len} bp   margins L={lm:.0f} R={rm:.0f}   "
                f"center_offset={ps.metrics['amplicon_center_offset_bp']:.1f} bp  GC={ps.metrics['amplicon_gc_pct']:.1f}%"
            )
            if background_used:
                def _sumline(prefix: str, label: str):
                    exact = "YES" if ps.metrics.get(f"{prefix}_has_exact_full_hit", 0) else "no"
                    mm_any = ps.metrics.get(f"{prefix}_offtarget_mm", 0)
                    mm_imp = ps.metrics.get(f"{prefix}_offtarget_mm_imp", 0)
                    src = ps.metrics.get(f"{prefix}_offtarget_from", "")
                    hit = ps.metrics.get(f"{prefix}_offtarget_hit", "")
                    print(f"  {label}: top={hit or '-'}  mm_any={int(mm_any)}  mm_imp={int(mm_imp)} ({src})  exact_full={exact}")

                _sumline("fp", "FP off-target (bg)")
                _sumline("rp", "RP off-target (bg)")
                _sumline("crrna", "crRNA off-target (bg)")

            if ps.flags:
                print(f"  Flags: {ps.flags}")

        # CSV (already sorted and ranked)
        write_pairs_csv(pairs_sorted, out_csv, tname, background_used)

        print(
            "\n[OK] Features + conservation + scoring complete. "
            + ("Included BLAST specificity vs background." if background_used
               else "No/failed background; used within-target mismatch profile for specificity.")
            + f" Ranked results written to: {out_csv}"
        )

    except subprocess.CalledProcessError as e:
        print("\nERROR running an external tool (makeblastdb/blastn or RNAfold).", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)
