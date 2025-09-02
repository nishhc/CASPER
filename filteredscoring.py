#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
score_pairs.py

End-to-end pipeline for primer/crRNA set generation + staged filtering + scoring.

Changes vs prior:
- Prefilter target defaults to 1000 (diverse subset before heavy features).
- Strong tie-resistance in scoring:
  * Two learned composites (entropy-weights and PCA-weights),
  * Rank/Borda composite,
  * Deterministic micro tiebreaker from sequence/coordinates.
- Background FASTA remains OPTIONAL; BLAST runs if provided & tools found,
  else falls back to on-target mismatch scanning.

USAGE
-----
# One-shot (generate + score), no background
python3 score_pairs.py gen-score --target-fasta target.fa --cas-type cas12a

# One-shot with background (uses BLAST if available)
python3 score_pairs.py gen-score --target-fasta target.fa --cas-type cas12a \
  --background-fasta genome.fa --threads 8

# Generate only
python3 score_pairs.py generate --target-fasta target.fa --cas-type cas12a

# Score an existing features CSV
python3 score_pairs.py score --input-csv pairs_features.csv --output-csv pairs_scored.csv
"""

from __future__ import annotations

import argparse, os, sys, re, csv, tempfile, subprocess, math, json, random, hashlib
from dataclasses import dataclass, asdict
from typing import Tuple, Optional, List, Dict, Union

# ============================================================
# Small utils
# ============================================================

def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")
    return seq.translate(table)[::-1]

def gc_pct(seq: str) -> float:
    s = seq.upper()
    g = s.count("G"); c = s.count("C")
    atgc = sum(s.count(x) for x in "ACGT")
    return 100.0 * (g + c) / max(1, atgc)

def max_homopolymer_run(seq: str) -> int:
    if not seq: return 0
    run = best = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            run += 1; best = max(best, run)
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
    out = []; i = 0; L = len(needle)
    if L == 0: return out
    while True:
        j = hay.find(needle, i)
        if j == -1: break
        out.append(j); i = j + 1  # allow overlaps
    return out

def hamming(a: str, b: str) -> int:
    assert len(a) == len(b)
    return sum(1 for x,y in zip(a,b) if x != y)

def count_approx_matches(hay: str, pattern: str, max_mm: int) -> int:
    L = len(pattern); n = 0
    for i in range(0, len(hay) - L + 1):
        if hamming(hay[i:i+L], pattern) <= max_mm:
            n += 1
    return n

def which(cmd: str) -> Optional[str]:
    from shutil import which as _which
    return _which(cmd)

# ============================================================
# I/O + parsing
# ============================================================

def load_first_fasta_seq(path: str) -> Tuple[str, str]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA not found: {path}")
    name = None; seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].strip() or "target"
                else:
                    break
            else:
                seq_chunks.append(line)
    if name is None: raise ValueError(f"No FASTA header in: {path}")
    seq = "".join(seq_chunks).upper().replace("U","T")
    if not seq: raise ValueError(f"No sequence found under first header in: {path}")
    return name, seq

def parse_range_2ints(text: str, name: str) -> Tuple[int, int]:
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", text)
    if not m: raise argparse.ArgumentTypeError(f"{name} must be like '30-35'")
    a, b = int(m.group(1)), int(m.group(2))
    if a < 1 or b < 1 or a > b: raise argparse.ArgumentTypeError(f"{name} invalid: {a}-{b}")
    return (a, b)

def parse_range_2floats(text: str, name: str, lo_ok=0.0, hi_ok=100.0) -> Tuple[float, float]:
    m = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*-\s*([0-9]*\.?[0-9]+)\s*$", text)
    if not m: raise argparse.ArgumentTypeError(f"{name} must be like '35-55'")
    a, b = float(m.group(1)), float(m.group(2))
    if a < lo_ok or b > hi_ok or a > b: raise argparse.ArgumentTypeError(f"{name} invalid: {a}-{b}")
    return (a, b)

def validate_cas(t: str) -> str:
    t2 = t.strip().lower()
    if t2 not in {"cas12a", "cas9", "cas13"}:
        raise argparse.ArgumentTypeError("cas-type must be one of: cas12a, cas9, cas13")
    return t2

def default_crrna_len_for(cas_type: str) -> int:
    return {"cas12a": 23, "cas9": 20, "cas13": 28}[cas_type]

# ============================================================
# Config
# ============================================================

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

# ============================================================
# Guides
# ============================================================

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
    t = target.upper(); L = len(t); out: List[Guide] = []
    # + strand PAM TTTV
    for i in range(0, L - 4 - guide_len + 1):
        pam = t[i:i+4]
        if pam.startswith("TTT") and pam[3] in "ACGT":
            prot_start = i + 4
            prot = t[prot_start: prot_start + guide_len]
            if len(prot) == guide_len:
                out.append(Guide("cas12a", "+", pam, i, prot, prot_start, gc_pct(prot), False))
    # - strand PAM AAAN on +
    for i in range(0, L - 4):
        pam_plus = t[i:i+4]
        if pam_plus.startswith("AAA") and pam_plus[3] in "ACGT":
            prot_end_plus = i
            prot_start_plus = prot_end_plus - guide_len
            if prot_start_plus >= 0:
                prot_plus = t[prot_start_plus: prot_end_plus]
                prot = revcomp(prot_plus)
                out.append(Guide("cas12a", "-", revcomp(pam_plus), i, prot, prot_start_plus, gc_pct(prot), False))
    return out

def find_crrna_sites_cas9(target: str, guide_len: int) -> List[Guide]:
    t = target.upper(); L = len(t); out: List[Guide] = []
    # + strand NGG
    for i in range(0, L - 2):
        pam = t[i:i+3]
        if len(pam) == 3 and pam[1:] == "GG":
            prot_end = i; prot_start = prot_end - guide_len
            if prot_start >= 0:
                prot = t[prot_start: prot_end]
                out.append(Guide("cas9", "+", pam, i, prot, prot_start, gc_pct(prot), False))
    # - strand CCN
    for i in range(0, L - 2):
        pam_plus = t[i:i+3]
        if len(pam_plus) == 3 and pam_plus[:2] == "CC":
            prot_start_plus = i + 3; prot_end_plus = prot_start_plus + guide_len
            if prot_end_plus <= L:
                prot_plus = t[prot_start_plus: prot_end_plus]
                prot = revcomp(prot_plus)
                out.append(Guide("cas9", "-", revcomp(pam_plus), i, prot, prot_start_plus, gc_pct(prot), False))
    return out

def find_crrna_sites_cas13(target: str, guide_len: int) -> List[Guide]:
    t = target.upper(); L = len(t); out: List[Guide] = []
    for i in range(0, L - guide_len + 1):
        prot = t[i:i+guide_len]
        out.append(Guide("cas13", "+", "", i, prot, i, gc_pct(prot), False))
    return out

def discover_crrnas(target_seq: str, cfg: Config) -> List[Guide]:
    cas = cfg.cas_type
    if cas == "cas12a": raw = find_crrna_sites_cas12a(target_seq, cfg.crrna_len)
    elif cas == "cas9": raw = find_crrna_sites_cas9(target_seq, cfg.crrna_len)
    else: raw = find_crrna_sites_cas13(target_seq, cfg.crrna_len)

    lo, hi = cfg.gc_crrna_range
    filtered: List[Guide] = []
    for g in raw:
        hp = has_homopolymer(g.protospacer_seq, cfg.max_homopolymer_crrna)
        g.has_homopolymer = hp
        if hp: continue
        if not (lo <= g.gc_pct <= hi): continue
        filtered.append(g)

    print(f"[CRRNA] Found {len(raw)} raw {cas} sites; {len(filtered)} passed (GC {lo}–{hi}%, homopolymer<{cfg.max_homopolymer_crrna}).")
    for g in filtered[:min(10,len(filtered))]:
        pam_info = f"PAM={g.pam_seq} at {g.pam_start}" if g.pam_seq else "PAM=NA"
        print(f"  {g.cas_type} {g.strand}  {pam_info}  prot@{g.protospacer_start}  GC={g.gc_pct:.1f}%  seq={g.protospacer_seq}")
    return filtered

# ============================================================
# Primers
# ============================================================

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
    t = target.upper(); out: List[Primer] = []
    lo_gc, hi_gc = cfg.gc_primer_range
    for L in range(cfg.primer_len_range[0], cfg.primer_len_range[1]+1):
        for i in range(0, len(t) - L + 1):
            s = t[i:i+L]; gc = gc_pct(s)
            if not (lo_gc <= gc <= hi_gc): continue
            if has_homopolymer(s, cfg.max_homopolymer): continue
            clamp = gentle_gc_clamp_ok(s)
            if clamp != 1: continue
            out.append(Primer(i, L, s, '+', gc, False, clamp))
    return out

def candidate_primers_minus_strand(target: str, cfg: Config) -> List[Primer]:
    t = target.upper(); out: List[Primer] = []
    lo_gc, hi_gc = cfg.gc_primer_range
    for L in range(cfg.primer_len_range[0], cfg.primer_len_range[1]+1):
        for i in range(0, len(t) - L + 1):
            win = t[i:i+L]; s = revcomp(win); gc = gc_pct(s)
            if not (lo_gc <= gc <= hi_gc): continue
            if has_homopolymer(s, cfg.max_homopolymer): continue
            clamp = gentle_gc_clamp_ok(s)
            if clamp != 1: continue
            out.append(Primer(i, L, s, '-', gc, False, clamp))
    return out

# ============================================================
# Pairing
# ============================================================

@dataclass
class PairedSet:
    guide_index: int
    guide: Guide
    fp: Primer
    rp: Primer
    amplicon_len: int
    metrics: Dict[str, float]
    flags: Dict[str, str]

def pair_primers(target: str, guides: List[Guide], fps: List[Primer], rps: List[Primer], cfg: Config) -> List[PairedSet]:
    pairs: List[PairedSet] = []
    min_amp, max_amp = cfg.amplicon_len_range; margin = cfg.guide_margin_bp
    for gi, g in enumerate(guides):
        left_ok = [fp for fp in fps if (fp.start + fp.length) <= (g.protospacer_start - margin)]
        guide_end = g.protospacer_start + cfg.crrna_len
        right_ok = [rp for rp in rps if rp.start >= (guide_end + margin)]
        for fp in left_ok:
            amplicon_start = fp.start
            for rp in right_ok:
                amplicon_end = rp.start + rp.length
                amplicon = amplicon_end - amplicon_start
                if not (min_amp <= amplicon <= max_amp): continue
                left_gap = g.protospacer_start - amplicon_start
                right_gap = amplicon_end - (g.protospacer_start + cfg.crrna_len)
                if left_gap < margin or right_gap < margin: continue
                pairs.append(PairedSet(gi, g, fp, rp, amplicon, {}, {}))
    return pairs

# ============================================================
# Cheap features (prefilter) and helper filters
# ============================================================

def annotate_light_features(ps: PairedSet, cfg: Config, target: str):
    """Cheap metrics: geometry, ΔTm (Wallace), amplicon/seed/primer GC, 3'-k uniqueness, on-target amplicon uniqueness."""
    amp_start = ps.fp.start
    amp_end = ps.rp.start + ps.rp.length
    amp_mid = amp_start + ps.amplicon_len / 2.0

    g_start = ps.guide.protospacer_start
    g_end = g_start + cfg.crrna_len
    g_mid = (g_start + g_end) / 2.0

    left_margin = g_start - amp_start
    right_margin = amp_end - g_end
    centered_abs_bp = abs(g_mid - amp_mid)

    seed_seq = ps.guide.protospacer_seq[-8:] if ps.guide.cas_type == "cas9" else ps.guide.protospacer_seq[:8]
    seed_gc = gc_pct(seed_seq)
    seed_run = max_homopolymer_run(seed_seq)

    amp_seq = target[amp_start:amp_end]
    amp_gc = gc_pct(amp_seq)
    amp_run = max_homopolymer_run(amp_seq)

    tm_fp = wallace_tm(ps.fp.seq)
    tm_rp = wallace_tm(ps.rp.seq)
    delta_tm = abs(tm_fp - tm_rp)

    # 3' k-mer uniqueness (on target)
    k = min(cfg.on_target_k, ps.fp.length, ps.rp.length)
    rp_site_top = revcomp(ps.rp.seq)
    fp_k = ps.fp.seq[-k:]
    rp_k = rp_site_top[-k:]
    tU = target.upper()
    fp_k_hits = len(find_all_indices(tU, fp_k))
    rp_k_hits = len(find_all_indices(tU, rp_k))

    # Amplicon uniqueness
    min_amp, max_amp = cfg.amplicon_len_range
    rp_len = len(ps.rp.seq)
    fp_starts = find_all_indices(tU, ps.fp.seq)
    rp_starts = find_all_indices(tU, rp_site_top)
    amplicon_hits = 0
    for i in fp_starts:
        for j in rp_starts:
            if j <= i: continue
            amp_len = (j + rp_len) - i
            if min_amp <= amp_len <= max_amp:
                amplicon_hits += 1

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

        "tm_fp_C": tm_fp,
        "tm_rp_C": tm_rp,
        "delta_tm_C": delta_tm,

        "fp_3p_k": float(k),
        "fp_3p_k_hits_target": float(fp_k_hits),
        "rp_3p_k_hits_target": float(rp_k_hits),

        "amplicon_hits_on_target": float(amplicon_hits),
    })

def passes_hard_light_filters(ps: PairedSet, cfg: Config) -> bool:
    if ps.metrics.get("left_margin_bp", 0) < cfg.guide_margin_bp: return False
    if ps.metrics.get("right_margin_bp", 0) < cfg.guide_margin_bp: return False
    if int(ps.metrics.get("amplicon_hits_on_target", 0)) != 1: return False
    if int(ps.metrics.get("fp_3p_k_hits_target", 2)) != 1: return False
    if int(ps.metrics.get("rp_3p_k_hits_target", 2)) != 1: return False
    if ps.metrics.get("delta_tm_C", 999) > 15.0: return False
    if int(ps.metrics.get("seed_has_homopolymer", 0)) == 1: return False
    return True

def light_score(ps: PairedSet, cfg: Config) -> float:
    """Fast 0..1 score for prefilter ranking; more continuous to reduce ties."""
    def cap01(x: float): return max(0.0, min(1.0, float(x)))
    def tri(x: float, lo: float, mid: float, hi: float):
        if x <= lo or x >= hi: return 0.0
        return (x - lo) / (mid - lo) if x < mid else (hi - x) / (hi - mid)

    lm = ps.metrics["left_margin_bp"]; rm = ps.metrics["right_margin_bp"]
    center_off = ps.metrics["amplicon_center_offset_bp"]
    geom = 0.4*cap01(lm / max(1.0, cfg.guide_margin_bp)) + \
           0.4*cap01(rm / max(1.0, cfg.guide_margin_bp)) + \
           0.2*(1.0/(1.0 + center_off/20.0))

    lo_amp, hi_amp = cfg.amplicon_len_range; mid_amp = (lo_amp+hi_amp)/2.0
    amp_len_part = tri(ps.amplicon_len, lo_amp, mid_amp, hi_amp)

    amp_gc = ps.metrics["amplicon_gc_pct"]; comp_amp = tri(amp_gc, 35.0, 50.0, 65.0)
    lo_gc_p, hi_gc_p = cfg.gc_primer_range; mid_gc_p = (lo_gc_p+hi_gc_p)/2.0
    comp_primer = 0.5*tri(ps.fp.gc_pct, lo_gc_p, mid_gc_p, hi_gc_p) + \
                  0.5*tri(ps.rp.gc_pct, lo_gc_p, mid_gc_p, hi_gc_p)

    seed_gc = ps.metrics["seed_gc_pct"]
    if ps.guide.cas_type == "cas9": sg_lo, sg_hi = (40.0, 60.0)
    else: sg_lo, sg_hi = (35.0, 55.0)
    seed_part = tri(seed_gc, sg_lo, (sg_lo+sg_hi)/2.0, sg_hi)

    dtm = ps.metrics["delta_tm_C"]; dtm_part = cap01(1.0 - dtm/15.0)

    # slight smooth preference for primer length balance (reduces ties)
    Ldiff = abs(ps.fp.length - ps.rp.length)
    len_bal = 1.0/(1.0 + Ldiff/2.0)

    s = (0.28*geom + 0.20*amp_len_part + 0.18*comp_amp +
         0.14*comp_primer + 0.12*seed_part + 0.06*dtm_part + 0.02*len_bal)
    return cap01(s)

def similar_key(ps: PairedSet, bin_bp: int) -> Tuple[int, int, int]:
    fp_bin = ps.fp.start // max(1, bin_bp)
    rp_bin = ps.rp.start // max(1, bin_bp)
    return (ps.guide_index, fp_bin, rp_bin)

# ============================================================
# Heavy features (after down-selection)
# ============================================================

def max_complementarity_percent(a: str, b: str) -> float:
    comp = str.maketrans("ACGT", "TGCA")
    aU = a.upper()
    bC = b.upper().translate(comp)[::-1]
    la, lb = len(aU), len(bC)
    if la == 0 or lb == 0: return 0.0
    max_match = 0
    for shift in range(-(lb-1), la):
        matches = 0; overlap = 0
        for i in range(lb):
            j = i + shift
            if 0 <= j < la:
                overlap += 1
                if bC[i] == aU[j]: matches += 1
        if overlap > 0: max_match = max(max_match, matches)
    denom = min(la, lb)
    return 100.0 * (max_match / max(1, denom))

def three_prime_run(a: str, b: str) -> int:
    aU = a.upper(); brc = revcomp(b.upper())
    m = min(len(aU), len(brc)); run = 0
    for k in range(1, m+1):
        if aU[-k] == brc[-k]: run += 1
        else: break
    return run

def three_prime_self_hairpin_run(seq: str) -> int:
    s = seq.upper()
    for k in range(min(12, len(s)), 0, -1):
        tail = s[-k:]; tail_rc = revcomp(tail)
        if tail_rc in s[:-k]: return k
    return 0

def has_gquad_seed_3p(seq: str, window: int = 12) -> bool:
    tail = seq[-window:].upper() if len(seq) >= window else seq.upper()
    return "GGGG" in tail

def annotate_structure_and_tm(ps: PairedSet, cfg: Config):
    tm_fp = wallace_tm(ps.fp.seq); tm_rp = wallace_tm(ps.rp.seq); delta_tm = abs(tm_fp - tm_rp)
    fp_self = max_complementarity_percent(ps.fp.seq, ps.fp.seq)
    rp_self = max_complementarity_percent(ps.rp.seq, ps.rp.seq)
    fp_rp_cross = max_complementarity_percent(ps.fp.seq, ps.rp.seq)
    fp_3p_self = three_prime_self_hairpin_run(ps.fp.seq)
    rp_3p_self = three_prime_self_hairpin_run(ps.rp.seq)
    fp_rp_3p = three_prime_run(ps.fp.seq, ps.rp.seq)
    fp_gquad = 1.0 if has_gquad_seed_3p(ps.fp.seq) else 0.0
    rp_gquad = 1.0 if has_gquad_seed_3p(ps.rp.seq) else 0.0

    if fp_self > cfg.self_dimer_thresh_pct: ps.flags["fp_self_dimer"] = f"{fp_self:.1f}%"
    if rp_self > cfg.self_dimer_thresh_pct: ps.flags["rp_self_dimer"] = f"{rp_self:.1f}%"
    if fp_rp_cross > cfg.cross_dimer_thresh_pct: ps.flags["fp_rp_cross_dimer"] = f"{fp_rp_cross:.1f}%"
    if delta_tm > cfg.tm_delta_warn: ps.flags["delta_tm_high"] = f"{delta_tm:.1f}C"
    if fp_rp_3p >= 4: ps.flags["fp_rp_3p_run_ge4"] = str(fp_rp_3p)
    if fp_3p_self >= 4: ps.flags["fp_3p_self_ge4"] = str(fp_3p_self)
    if rp_3p_self >= 4: ps.flags["rp_3p_self_ge4"] = str(rp_3p_self)

    ps.metrics.update({
        "tm_fp_C": tm_fp, "tm_rp_C": tm_rp, "delta_tm_C": delta_tm,
        "fp_self_dimer_pct": fp_self, "rp_self_dimer_pct": rp_self,
        "fp_rp_cross_dimer_pct": fp_rp_cross,
        "fp_3p_self_run": float(fp_3p_self), "rp_3p_self_run": float(rp_3p_self),
        "fp_rp_3p_cross_run": float(fp_rp_3p),
        "fp_gquad_3p": fp_gquad, "rp_gquad_3p": rp_gquad,
        "fp_max_hpoly_run": float(max_homopolymer_run(ps.fp.seq)),
        "rp_max_hpoly_run": float(max_homopolymer_run(ps.rp.seq)),
        "fp_tail_gc": float(ps.fp.tail_gc_count), "rp_tail_gc": float(ps.rp.tail_gc_count),
    })

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

# ============================================================
# ViennaRNA (optional; cached per guide)
# ============================================================

_MFE_CACHE: Dict[str, Tuple[float, str, float]] = {}

def rnafold_cli_fold(seq_rna: str) -> Tuple[Optional[float], Optional[str]]:
    if which("RNAfold") is None: return (None, None)
    try:
        p = subprocess.run(["RNAfold","--noPS"], input=(seq_rna+"\n").encode(),
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        out = p.stdout.decode().strip().splitlines()
        if len(out) >= 2:
            m = re.search(r"([().]+)\s+\(([-0-9.]+)\)", out[1])
            if m: return (float(m.group(2)), str(m.group(1)))
    except Exception:
        pass
    return (None, None)

def viennarna_mfe_fold(seq_rna: str) -> Tuple[Optional[float], Optional[str]]:
    try:
        import RNA
    except Exception:
        return (None, None)
    try:
        fc = RNA.fold_compound(seq_rna); dot, mfe = fc.mfe()
        return (float(mfe), str(dot))
    except Exception:
        return (None, None)

def rna_secondary_structure(seq_dna: str) -> Tuple[Optional[float], Optional[str]]:
    seq_rna = seq_dna.replace("T","U")
    mfe, dot = viennarna_mfe_fold(seq_rna)
    if mfe is not None and dot is not None and len(dot) == len(seq_dna): return (mfe, dot)
    mfe, dot = rnafold_cli_fold(seq_rna)
    if mfe is not None and dot is not None and len(dot) == len(seq_dna): return (mfe, dot)
    return (None, None)

def annotate_guide_secondary_structure(ps: PairedSet, cfg: Config, skip_vienna: bool=False):
    g = ps.guide.protospacer_seq
    if skip_vienna:
        ps.metrics.update({"guide_mfe_kcal": 0.0, "guide_dotbracket": "", "guide_seed_unpaired_frac": 0.5})
        return
    if g in _MFE_CACHE:
        mfe, dot, frac = _MFE_CACHE[g]
        ps.metrics.update({"guide_mfe_kcal": mfe, "guide_dotbracket": dot, "guide_seed_unpaired_frac": frac})
        return
    mfe, dot = rna_secondary_structure(g)
    seed_slice = slice(-8, None) if ps.guide.cas_type == "cas9" else slice(0, 8)
    if dot and len(dot) == len(g):
        seed = dot[seed_slice]; frac = seed.count(".") / max(1, len(seed))
        mfe_out = float(mfe) if mfe is not None else 0.0
        _MFE_CACHE[g] = (mfe_out, dot, float(frac))
        ps.metrics.update({"guide_mfe_kcal": mfe_out, "guide_dotbracket": dot, "guide_seed_unpaired_frac": float(frac)})
    else:
        _MFE_CACHE[g] = (0.0, "", 0.5)
        ps.metrics.update({"guide_mfe_kcal": 0.0, "guide_dotbracket": "", "guide_seed_unpaired_frac": 0.5})

# ============================================================
# Specificity: BLAST vs on-target fallback
# ============================================================

Hit = Dict[str, Union[float, int, str]]

class BlastCache:
    def __init__(self): self.cache: Dict[str, List[Hit]] = {}
    def get(self, seq: str): return self.cache.get(seq)
    def set(self, seq: str, val: List[Hit]): self.cache[seq] = val

def ensure_blast_db(background_fasta: str, workdir: str) -> str:
    if which("makeblastdb") is None:
        raise RuntimeError("makeblastdb not found in PATH.")
    db_prefix = os.path.join(workdir, "background_db")
    cmd = ["makeblastdb","-in",background_fasta,"-dbtype","nucl","-parse_seqids","-out",db_prefix]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return db_prefix

def blast_collect_hits(query_seq: str, db_prefix: str, cfg: Config, workdir: str) -> List[Hit]:
    cachedir = os.path.join(workdir, "q"); os.makedirs(cachedir, exist_ok=True)
    qfa = os.path.join(cachedir, "q.fa"); out_tsv = os.path.join(cachedir, "hits.tsv")
    with open(qfa, "w") as f:
        f.write(">q\n"); f.write(query_seq + "\n")
    cmd = [
        "blastn","-query",qfa,"-db",db_prefix,
        "-word_size",str(cfg.blast_word_size),
        "-reward","1","-penalty","-2","-gapopen","5","-gapextend","2",
        "-evalue",str(cfg.blast_evalue),
        "-strand","both","-task",cfg.blast_task,
        "-dust","no","-soft_masking","false",
        "-outfmt","6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore sstrand",
        "-out",out_tsv
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    hits: List[Hit] = []
    qlen = int(len(query_seq))
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) < 12: continue
                pident = float(p[2]); aln_len = int(p[3]); mismatch = int(p[4])
                qstart = int(p[5]); qend = int(p[6]); sstart = int(p[7]); send = int(p[8])
                sseqid = p[1]; sstrand = p[11]
                score = (pident/100.0) * (aln_len / max(1.0, float(qlen)))
                is_exact_full = (mismatch == 0 and aln_len == qlen and ((qstart == 1 and qend == qlen) or (qstart == qlen and qend == 1)))
                hits.append({
                    "sseqid": sseqid, "pident": pident, "aln_len": aln_len, "mismatch": mismatch,
                    "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
                    "sstrand": sstrand, "score": score, "is_exact_full": 1 if is_exact_full else 0
                })
    return hits

def _top_any_and_imperfect(hits: List[Hit], qlen: int):
    if not hits: return None, None, 0, ""
    hits_sorted = sorted(hits, key=lambda h: h["score"], reverse=True)
    top_any = hits_sorted[0]
    n_exact = sum(1 for h in hits_sorted if int(h["is_exact_full"]) == 1)
    example_exact_id = ""
    for h in hits_sorted:
        if int(h["is_exact_full"]) == 1:
            example_exact_id = str(h["sseqid"]); break
    imperfect = None
    for h in hits_sorted:
        if int(h["mismatch"]) > 0 or int(h["aln_len"]) != qlen:
            imperfect = h; break
    return top_any, imperfect, n_exact, example_exact_id

def annotate_specificity_with_blast(pairs: List[PairedSet], cfg: Config, background_fasta: str):
    if which("blastn") is None: raise RuntimeError("blastn not found in PATH.")
    with tempfile.TemporaryDirectory(prefix="primedrpa_blast_") as td:
        db_prefix = ensure_blast_db(background_fasta, td)
        cache = BlastCache()
        def get_hits(seq: str):
            cached = cache.get(seq)
            if cached is not None: return cached
            h = blast_collect_hits(seq, db_prefix, cfg, td)
            cache.set(seq, h); return h
        for ps in pairs:
            fp_hits = get_hits(ps.fp.seq); rp_hits = get_hits(ps.rp.seq); g_hits = get_hits(ps.guide.protospacer_seq)
            qlen_fp = len(ps.fp.seq); qlen_rp = len(ps.rp.seq); qlen_g = len(ps.guide.protospacer_seq)
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

def annotate_specificity_on_target_fallback(ps: PairedSet, target: str):
    """Populate off-target metrics by scanning the target itself, excluding intended window."""
    t = target.upper()
    def best_other_window(pattern: str, intended_start: int) -> Dict[str, Union[float,int,str]]:
        L = len(pattern); best_ident = -1.0; best_mm = L; best_i = -1
        for i in range(0, len(t) - L + 1):
            if i == intended_start: continue
            mm = hamming(t[i:i+L], pattern); ident = (L - mm) / float(L)
            if ident > best_ident:
                best_ident = ident; best_mm = mm; best_i = i
        exact_positions = [i for i in range(0, len(t) - L + 1) if t[i:i+L] == pattern]
        exact_other = max(0, len(exact_positions) - (1 if intended_start in exact_positions else 0))
        return {
            "score": float(max(0.0, best_ident)),
            "mismatch": int(best_mm),
            "aln_len": int(L),
            "hit": f"target:{best_i}" if best_i >= 0 else "",
            "has_exact_other": 1 if exact_other > 0 else 0,
            "exact_other_count": int(exact_other)
        }
    fp_res = best_other_window(ps.fp.seq, ps.fp.start)
    rp_site_top = revcomp(ps.rp.seq)
    rp_res = best_other_window(rp_site_top, ps.rp.start)
    g = ps.guide
    g_res = best_other_window(g.protospacer_seq, g.protospacer_start)
    def save(prefix: str, res: Dict[str,Union[float,int,str]]):
        if res["hit"] == "":
            ps.metrics[f"{prefix}_offtarget_score"] = 0.0
            ps.metrics[f"{prefix}_offtarget_hit"] = ""
            ps.metrics[f"{prefix}_offtarget_mm"] = int(res["aln_len"])
            ps.metrics[f"{prefix}_offtarget_aln"] = 0
            ps.metrics[f"{prefix}_offtarget_score_imp"] = 0.0
            ps.metrics[f"{prefix}_offtarget_mm_imp"] = int(res["aln_len"])
            ps.metrics[f"{prefix}_offtarget_from"] = "none"
        else:
            ps.metrics[f"{prefix}_offtarget_score"] = float(res["score"])
            ps.metrics[f"{prefix}_offtarget_hit"] = str(res["hit"])
            ps.metrics[f"{prefix}_offtarget_mm"] = int(res["mismatch"])
            ps.metrics[f"{prefix}_offtarget_aln"] = int(res["aln_len"])
            ps.metrics[f"{prefix}_offtarget_score_imp"] = float(res["score"])
            ps.metrics[f"{prefix}_offtarget_mm_imp"] = int(res["mismatch"])
            ps.metrics[f"{prefix}_offtarget_from"] = "imperfect" if int(res["mismatch"])>0 else "any"
        ps.metrics[f"{prefix}_has_exact_full_hit"] = int(res["has_exact_other"])
        ps.metrics[f"{prefix}_exact_full_hit_count"] = int(res["exact_other_count"])
        ps.metrics[f"{prefix}_exact_full_example_sseqid"] = "target" if int(res["exact_other_count"])>0 else ""
    save("fp", fp_res); save("rp", rp_res); save("crrna", g_res)

# ============================================================
# CSV writer
# ============================================================

def write_pairs_csv(pairs: List[PairedSet], out_path: str, target_name: str):
    fields = [
        "target", "cas_type", "guide_idx", "guide_strand",
        "guide_pam", "guide_pam_start", "guide_start", "guide_len",
        "guide_gc_pct", "guide_seq",
        "fp_start", "fp_len", "fp_gc_pct", "fp_seq", "fp_tail_gc", "fp_max_hpoly_run",
        "rp_start", "rp_len", "rp_gc_pct", "rp_seq", "rp_tail_gc", "rp_max_hpoly_run",
        "amplicon_len",
        "left_margin_bp", "right_margin_bp", "amplicon_center_offset_bp",
        "tm_fp_C", "tm_rp_C", "delta_tm_C",
        "fp_self_dimer_pct", "rp_self_dimer_pct", "fp_rp_cross_dimer_pct",
        "fp_3p_self_run", "rp_3p_self_run", "fp_rp_3p_cross_run",
        "fp_gquad_3p", "rp_gquad_3p",
        "seed_seq", "seed_len", "seed_gc_pct", "seed_max_run", "seed_has_homopolymer",
        "amplicon_gc_pct", "amplicon_max_run",
        "overlap_protospacer", "overlap_pam",
        "fp_guide_cross_pct", "rp_guide_cross_pct", "fp_guide_3p_run", "rp_guide_3p_run",
        "fp_full_hits_target_both", "rp_full_hits_target_both",
        "fp_full_hits_mm1_both", "rp_full_hits_mm1_both",
        "fp_3p_k", "fp_3p_k_hits_target", "rp_3p_k_hits_target",
        "fp_3p_k_hits_mm1", "rp_3p_k_hits_mm1",
        "amplicon_hits_on_target",
        "fp_offtarget_score", "fp_offtarget_hit", "fp_offtarget_mm", "fp_offtarget_aln",
        "fp_offtarget_score_imp", "fp_offtarget_mm_imp", "fp_offtarget_from",
        "fp_has_exact_full_hit", "fp_exact_full_hit_count", "fp_exact_full_example_sseqid",
        "rp_offtarget_score", "rp_offtarget_hit", "rp_offtarget_mm", "rp_offtarget_aln",
        "rp_offtarget_score_imp", "rp_offtarget_mm_imp", "rp_offtarget_from",
        "rp_has_exact_full_hit", "rp_exact_full_hit_count", "rp_exact_full_example_sseqid",
        "crrna_offtarget_score", "crrna_offtarget_hit", "crrna_offtarget_mm", "crrna_offtarget_aln",
        "crrna_offtarget_score_imp", "crrna_offtarget_mm_imp", "crrna_offtarget_from",
        "crrna_has_exact_full_hit", "crrna_exact_full_hit_count", "crrna_exact_full_example_sseqid",
        "guide_mfe_kcal", "guide_seed_unpaired_frac", "guide_dotbracket",
        "flags"
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for ps in pairs:
            g = ps.guide
            row = {
                "target": target_name, "cas_type": g.cas_type,
                "guide_idx": ps.guide_index, "guide_strand": g.strand,
                "guide_pam": g.pam_seq, "guide_pam_start": g.pam_start,
                "guide_start": g.protospacer_start, "guide_len": len(g.protospacer_seq),
                "guide_gc_pct": f"{g.gc_pct:.1f}", "guide_seq": g.protospacer_seq,
                "fp_start": ps.fp.start, "fp_len": ps.fp.length, "fp_gc_pct": f"{ps.fp.gc_pct:.1f}",
                "fp_seq": ps.fp.seq, "fp_tail_gc": int(ps.metrics.get("fp_tail_gc", 0)),
                "fp_max_hpoly_run": int(ps.metrics.get("fp_max_hpoly_run", 0)),
                "rp_start": ps.rp.start, "rp_len": ps.rp.length, "rp_gc_pct": f"{ps.rp.gc_pct:.1f}",
                "rp_seq": ps.rp.seq, "rp_tail_gc": int(ps.metrics.get("rp_tail_gc", 0)),
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
                "seed_seq": ps.metrics.get("seed_seq", ""), "seed_len": int(ps.metrics.get("seed_len", 0)),
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
                "fp_full_hits_mm1_both": int(ps.metrics.get('fp_full_hits_mm1_both', 0)),
                "rp_full_hits_mm1_both": int(ps.metrics.get('rp_full_hits_mm1_both', 0)),
                "fp_3p_k": int(ps.metrics.get("fp_3p_k", 0)),
                "fp_3p_k_hits_target": int(ps.metrics.get('fp_3p_k_hits_target', 0)),
                "rp_3p_k_hits_target": int(ps.metrics.get('rp_3p_k_hits_target', 0)),
                "fp_3p_k_hits_mm1": int(ps.metrics.get('fp_3p_k_hits_mm1', 0)),
                "rp_3p_k_hits_mm1": int(ps.metrics.get('rp_3p_k_hits_mm1', 0)),
                "amplicon_hits_on_target": int(ps.metrics.get('amplicon_hits_on_target', 0)),
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
                "flags": ";".join([f"{k}={v}" for k, v in ps.flags.items()]) if ps.flags else ""
            }
            w.writerow(row)
    print(f"[WRITE] CSV → {out_path}  ({len(pairs)} rows)")

# ============================================================
# Generator main with staged prefilter & diversity
# ============================================================

def generator_main(args) -> str:
    target_name, target_seq = load_first_fasta_seq(args.target_fasta)
    bg_path = None
    if args.background_fasta:
        if os.path.exists(args.background_fasta):
            bg_path = os.path.abspath(args.background_fasta)
        else:
            print(f"[WARN] --background-fasta not found ({args.background_fasta}); proceeding without background.", file=sys.stderr)

    cas = args.cas_type
    crrna_len = args.crrna_len if args.crrna_len is not None else default_crrna_len_for(cas)

    cfg = Config(
        cas_type=cas, top_k=max(1, int(args.top_k)),
        primer_len_range=args.primer_len, amplicon_len_range=args.amplicon_len,
        gc_primer_range=args.gc_primer, max_homopolymer=max(2, int(args.max_homopolymer)),
        crrna_len=int(crrna_len), gc_crrna_range=args.gc_crrna, max_homopolymer_crrna=max(2, int(args.max_homopolymer_crrna)),
        cross_dimer_thresh_pct=float(args.cross_dimer_thresh), self_dimer_thresh_pct=float(args.self_dimer_thresh),
        tm_delta_warn=float(args.tm_delta_warn),
        blast_word_size=int(args.blast_word_size), blast_evalue=float(args.blast_evalue),
        blast_task=str(args.blast_task), threads=max(1, int(args.threads)),
        target_name=target_name, background_fasta=bg_path,
        guide_margin_bp=max(0, int(args.guide_margin)),
        on_target_k=max(6, int(args.on_target_k)),
        offtarget_max_allow=float(args.offtarget_max_allow),
    )

    if args.random_seed is not None:
        random.seed(int(args.random_seed))

    print("\n[CONFIG]")
    for k, v in asdict(cfg).items(): print(f"  {k}: {v}")
    print(f"\n[INPUT] target {target_name}: {len(target_seq)} bp")
    print(f"[INPUT] background FASTA: {cfg.background_fasta or '(none; on-target fallback)'}")
    print(f"[PREFILTER] target={args.prefilter_target}  similarity-bin={args.similarity_bin} bp  per-guide-cap-mult={args.per_guide_cap_mult}\n")

    # Guides / primers
    guides = discover_crrnas(target_seq, cfg)
    if not guides:
        print("\n[NOTE] No crRNA candidates passed filtering."); return ""

    fps = candidate_primers_plus_strand(target_seq, cfg)
    rps = candidate_primers_minus_strand(target_seq, cfg)
    print(f"\n[PRIMERS] {len(fps)} forward; {len(rps)} reverse candidates passed cheap filters "
          f"(GC {cfg.gc_primer_range[0]}–{cfg.gc_primer_range[1]}%, homopolymer<{cfg.max_homopolymer}, 3' clamp=1 GC).")
    if not fps or not rps:
        print("[NOTE] No primer candidates after cheap filters."); return ""

    # Pair
    pairs_all = pair_primers(target_seq, guides, fps, rps, cfg)
    print(f"[PAIRING] {len(pairs_all)} FP/RP pairs span guides with margin {cfg.guide_margin_bp} bp and "
          f"amplicon {cfg.amplicon_len_range[0]}–{cfg.amplicon_len_range[1]} bp.")
    if not pairs_all:
        print("[NOTE] No pairs satisfied geometry/amplicon constraints."); return ""

    # Stage 0: cheap features + hard light filters
    for ps in pairs_all:
        annotate_light_features(ps, cfg, target_seq)
    pairs_light = [ps for ps in pairs_all if passes_hard_light_filters(ps, cfg)]
    print(f"[PREFILTER] Hard light filters kept {len(pairs_light)} / {len(pairs_all)}.")
    if not pairs_light:
        print("[NOTE] Nothing left after prefilter."); return ""

    # Rank by light score
    scored = [(light_score(ps, cfg), idx, ps) for idx, ps in enumerate(pairs_light)]
    scored.sort(key=lambda z: z[0], reverse=True)

    # Greedy diversity selection
    target_n = max(1, int(args.prefilter_target))
    gcount = len({ps.guide_index for _,__,ps in scored})
    base_quota = max(1, target_n // max(1, gcount))
    max_per_guide = int(base_quota * float(args.per_guide_cap_mult)) + 2

    used_bins: set = set()
    taken: List[PairedSet] = []
    per_guide: Dict[int, int] = {}

    for s, _, ps in scored:
        if len(taken) >= target_n: break
        key = similar_key(ps, args.similarity_bin)
        if key in used_bins: continue
        if per_guide.get(ps.guide_index, 0) >= max_per_guide: continue
        used_bins.add(key)
        per_guide[ps.guide_index] = per_guide.get(ps.guide_index, 0) + 1
        taken.append(ps)

    if len(taken) < target_n:
        for s, _, ps in scored:
            if len(taken) >= target_n: break
            key = similar_key(ps, args.similarity_bin)
            if key in used_bins: continue
            used_bins.add(key)
            per_guide[ps.guide_index] = per_guide.get(ps.guide_index, 0) + 1
            taken.append(ps)

    if len(taken) < target_n and args.similarity_bin > 1:
        loosen = max(1, args.similarity_bin - 1)
        used_bins = set(); retaken = []; per_guide = {}
        for s, _, ps in scored:
            if len(retaken) >= target_n: break
            key = similar_key(ps, loosen)
            if key in used_bins: continue
            used_bins.add(key)
            per_guide[ps.guide_index] = per_guide.get(ps.guide_index, 0) + 1
            retaken.append(ps)
        if len(retaken) > len(taken):
            taken = retaken

    print(f"[PREFILTER] Diversity selection kept {len(taken)} (target {target_n}); "
          f"guides represented={len({ps.guide_index for ps in taken})}.")

    # Stage 1: heavy features on the kept set only
    for ps in taken:
        amp_start = ps.fp.start; amp_end = ps.rp.start + ps.rp.length
        g_start = ps.guide.protospacer_start; g_end = g_start + cfg.crrna_len
        pam_start = ps.guide.pam_start; pam_len = len(ps.guide.pam_seq) if ps.guide.pam_seq else 0
        pam_end = pam_start + pam_len
        fp_span = (ps.fp.start, ps.fp.start + ps.fp.length)
        rp_span = (ps.rp.start, ps.rp.start + ps.rp.length)
        guide_span = (g_start, g_end)
        def spans_overlap(a: Tuple[int,int], b: Tuple[int,int]) -> bool:
            return not (a[1] <= b[0] or b[1] <= a[0])
        overlap_protospacer = (spans_overlap(fp_span, guide_span) or spans_overlap(rp_span, guide_span))
        overlap_pam = (pam_len > 0) and (spans_overlap(fp_span, (pam_start, pam_end)) or spans_overlap(rp_span, (pam_start, pam_end)))
        ps.metrics.update({"overlap_protospacer": 1.0 if overlap_protospacer else 0.0,
                           "overlap_pam": 1.0 if overlap_pam else 0.0})
        if overlap_protospacer: ps.flags["primer_overlaps_protospacer"] = "1"
        if overlap_pam: ps.flags["primer_overlaps_pam"] = "1"

        annotate_structure_and_tm(ps, cfg)
        annotate_cross_with_guide(ps)
        annotate_guide_secondary_structure(ps, cfg, skip_vienna=bool(args.skip_vienna))

    # Specificity: BLAST if background+tools else fast fallback
    if cfg.background_fasta and not args.skip_blast:
        try:
            print("[BLAST] Using background DB for FP, RP, and crRNA (per-component, mismatch-aware)...")
            annotate_specificity_with_blast(taken, cfg, cfg.background_fasta)
        except Exception as e:
            print(f"[WARN] BLAST unavailable or failed ({e}). Falling back to on-target scan.", file=sys.stderr)
            for ps in taken:
                annotate_specificity_on_target_fallback(ps, target_seq)
    else:
        print("[SPEC] Background not used (or --skip-blast). Using on-target scan fallback for FP/RP/crRNA specificity.")
        for ps in taken:
            annotate_specificity_on_target_fallback(ps, target_seq)

    # Console preview
    show_n = min(cfg.top_k, len(taken))
    for ps in taken[:show_n]:
        g = ps.guide
        print(f"\n[SET] GuideIdx={ps.guide_index} {g.cas_type}{g.strand} "
              f"{'PAM='+g.pam_seq if g.pam_seq else 'PAM=NA'} prot@{g.protospacer_start} len={len(g.protospacer_seq)} GC={g.gc_pct:.1f}%")
        print(f"  FP  start={ps.fp.start} len={ps.fp.length} GC={ps.fp.gc_pct:.1f}%  Tm={ps.metrics['tm_fp_C']:.1f}°C  self%={ps.metrics.get('fp_self_dimer_pct',0):.1f}  seq={ps.fp.seq}")
        print(f"  RP  start={ps.rp.start} len={ps.rp.length} GC={ps.rp.gc_pct:.1f}%  Tm={ps.metrics['tm_rp_C']:.1f}°C  self%={ps.metrics.get('rp_self_dimer_pct',0):.1f}  seq={ps.rp.seq}")
        print(f"  Cross-dimer%={ps.metrics.get('fp_rp_cross_dimer_pct',0):.1f}   ΔTm={ps.metrics['delta_tm_C']:.1f}°C   3'cross_run={int(ps.metrics.get('fp_rp_3p_cross_run',0))}")
        lm = ps.metrics['left_margin_bp']; rm = ps.metrics['right_margin_bp']
        print(f"  Amplicon={ps.amplicon_len} bp   margins L={lm:.0f} R={rm:.0f}   center_offset={ps.metrics['amplicon_center_offset_bp']:.1f} bp  GC={ps.metrics['amplicon_gc_pct']:.1f}%")
        print(f"  Seed seq={ps.metrics['seed_seq']}  GC={ps.metrics['seed_gc_pct']:.1f}%  max_run={int(ps.metrics['seed_max_run'])}")
        def _sumline(prefix: str, label: str):
            exact = 'YES' if ps.metrics.get(f"{prefix}_has_exact_full_hit",0) else 'no'
            mm_any = ps.metrics.get(f"{prefix}_offtarget_mm",0)
            mm_imp = ps.metrics.get(f"{prefix}_offtarget_mm_imp",0)
            src = ps.metrics.get(f"{prefix}_offtarget_from","")
            hit = ps.metrics.get(f"{prefix}_offtarget_hit","")
            print(f"  {label}: top={hit or '-'}  mm_any={int(mm_any)}  mm_imp={int(mm_imp)} ({src})  extra_exact_elsewhere={exact}")
        _sumline("fp", "FP off-target"); _sumline("rp", "RP off-target"); _sumline("crrna", "crRNA off-target")
        if ps.flags: print(f"  Flags: {ps.flags}")

    # Write CSV
    write_pairs_csv(taken, args.write_csv, target_name)
    return args.write_csv

# ============================================================
# ---- Scoring (with strong tie resistance) ----
# ============================================================

try:
    import numpy as np
except Exception:
    np = None

def cap01(x: float) -> float: return max(0.0, min(1.0, float(x)))
def tri(x: float, lo: float, mid: float, hi: float):
    if x <= lo or x >= hi: return 0.0
    return (x - lo) / (mid - lo) if x < mid else (hi - x) / (hi - mid)
def inv_small(x: float, soft_ref: float) -> float: return cap01(1.0 - (x / max(1e-9, soft_ref)))
def inv_pct(x: float, soft_ref: float) -> float: return cap01(1.0 - (x / max(1e-9, soft_ref)))
def inv_pct_100(x: float) -> float: return cap01(1.0 - x / 100.0)
def norm_offtarget(score: float) -> float: return cap01(1.0 - score)
def mfe_soft_good(mfe: Optional[float]) -> float:
    if mfe is None: return 0.5
    if mfe >= -2.0: return 1.0
    if mfe <= -12.0: return 0.0
    return (mfe + 12.0) / 10.0

def ffloat(d: Dict[str,str], key: str, default: float=0.0) -> float:
    v = d.get(key, "")
    if v is None: return default
    v = str(v).strip()
    if v == "" or v.upper() in {"NA","NAN","NULL"}: return default
    try: return float(v)
    except Exception: return default

def fint(d: Dict[str,str], key: str, default: int=0) -> int:
    try: return int(round(ffloat(d, key, float(default))))
    except Exception: return default

class SParams:
    def __init__(self, gc_primer, amplicon_len, guide_margin_bp, tm_delta_warn, cross_dimer_thresh, self_dimer_thresh, seed_gc_cas9, seed_gc_other):
        self.gc_primer = gc_primer
        self.amplicon_len = amplicon_len
        self.guide_margin_bp = guide_margin_bp
        self.tm_delta_warn = tm_delta_warn
        self.cross_dimer_thresh = cross_dimer_thresh
        self.self_dimer_thresh = self_dimer_thresh
        self.seed_gc_cas9 = seed_gc_cas9
        self.seed_gc_other = seed_gc_other

def build_feature_scores_from_row(row: Dict[str,str], P: SParams) -> Tuple[List[str], List[float]]:
    names: List[str] = []; vals: List[float] = []
    # Specificity (invert similarity)
    names += ["spec_fp", "spec_rp", "spec_crrna"]
    vals  += [norm_offtarget(ffloat(row,"fp_offtarget_score",0.0)),
              norm_offtarget(ffloat(row,"rp_offtarget_score",0.0)),
              norm_offtarget(ffloat(row,"crrna_offtarget_score",0.0)),]
    # Exact off-target presence → goodness
    names += ["spec_no_exact_fp","spec_no_exact_rp","spec_no_exact_crrna"]
    vals  += [1.0 - cap01(ffloat(row,"fp_has_exact_full_hit",0.0)),
              1.0 - cap01(ffloat(row,"rp_has_exact_full_hit",0.0)),
              1.0 - cap01(ffloat(row,"crrna_has_exact_full_hit",0.0)),]
    # Dimerization / ΔTm
    names += ["cross_dimer","fp_self_dimer","rp_self_dimer","delta_tm"]
    vals  += [inv_pct(ffloat(row,"fp_rp_cross_dimer_pct",0.0), P.cross_dimer_thresh),
              inv_pct(ffloat(row,"fp_self_dimer_pct",0.0),     P.self_dimer_thresh),
              inv_pct(ffloat(row,"rp_self_dimer_pct",0.0),     P.self_dimer_thresh),
              inv_small(ffloat(row,"delta_tm_C",0.0),          P.tm_delta_warn)]
    # 3′ complementarity runs
    names += ["fp_rp_3p_cross_run","fp_3p_self_run","rp_3p_self_run"]
    vals  += [inv_small(ffloat(row,"fp_rp_3p_cross_run",0.0),4.0),
              inv_small(ffloat(row,"fp_3p_self_run",0.0),    4.0),
              inv_small(ffloat(row,"rp_3p_self_run",0.0),    4.0)]
    # 3' G-quad flags (invert)
    names += ["fp_no_gquad3p","rp_no_gquad3p"]
    vals  += [1.0 - cap01(ffloat(row,"fp_gquad_3p",0.0)),
              1.0 - cap01(ffloat(row,"rp_gquad_3p",0.0))]
    # Geometry
    names += ["left_margin","right_margin","amplicon_len"]
    lo_amp, hi_amp = P.amplicon_len; mid_amp = (lo_amp+hi_amp)/2.0
    vals  += [min(1.0, ffloat(row,"left_margin_bp",0.0)/max(1.0,P.guide_margin_bp)),
              min(1.0, ffloat(row,"right_margin_bp",0.0)/max(1.0,P.guide_margin_bp)),
              tri(ffloat(row,"amplicon_len",0.0), lo_amp, mid_amp, hi_amp)]
    # Uniqueness
    names += ["amplicon_unique"]; vals += [1.0 if fint(row,"amplicon_hits_on_target",1) == 1 else 0.0]
    # Composition
    lo_gc_p, hi_gc_p = P.gc_primer; mid_gc_p = (lo_gc_p+hi_gc_p)/2.0
    names += ["fp_gc","rp_gc"]
    vals  += [tri(ffloat(row,"fp_gc_pct",0.0), lo_gc_p, mid_gc_p, hi_gc_p),
              tri(ffloat(row,"rp_gc_pct",0.0), lo_gc_p, mid_gc_p, hi_gc_p)]
    # Seed GC (cas-specific)
    seed_gc = ffloat(row,"seed_gc_pct",50.0)
    cas = (row.get("cas_type","") or "").lower()
    lo_seed, hi_seed = (P.seed_gc_cas9 if cas == "cas9" else P.seed_gc_other)
    names += ["seed_gc"]; vals += [tri(seed_gc, lo_seed, (lo_seed+hi_seed)/2.0, hi_seed)]
    # Seed homopolymer
    names += ["seed_no_hpoly"]; vals += [1.0 - cap01(ffloat(row,"seed_has_homopolymer",0.0))]
    # Complexity
    names += ["amplicon_low_hpoly","fp_low_hpoly","rp_low_hpoly"]
    vals  += [inv_small(ffloat(row,"amplicon_max_run",1.0),8.0),
              inv_small(ffloat(row,"fp_max_hpoly_run",1.0),6.0),
              inv_small(ffloat(row,"rp_max_hpoly_run",1.0),6.0)]
    # Cross-hyb with guide
    names += ["fp_low_guide_cross","rp_low_guide_cross"]
    vals  += [inv_pct_100(ffloat(row,"fp_guide_cross_pct",0.0)),
              inv_pct_100(ffloat(row,"rp_guide_cross_pct",0.0))]
    # Guide accessibility
    names += ["guide_seed_access","guide_soft_mfe"]
    vals  += [cap01(ffloat(row,"guide_seed_unpaired_frac",0.5)),
              mfe_soft_good(ffloat(row,"guide_mfe_kcal",0.0))]
    vals = [cap01(v) for v in vals]
    return names, vals

def entropy_weights(X: "np.ndarray") -> "np.ndarray":
    eps = 1e-12; col_sum = X.sum(axis=0); col_sum[col_sum==0] = eps
    P = X / col_sum; n = max(2, X.shape[0]); k = 1.0 / math.log(n)
    H = -k * (P * np.log(P + eps)).sum(axis=0); D = 1.0 - H; D[D<0] = 0
    s = D.sum(); return D / (s if s > eps else 1.0)

def pca_weights(X: "np.ndarray") -> "np.ndarray":
    X0 = X - X.mean(axis=0, keepdims=True)
    s = X0.std(axis=0, ddof=1, keepdims=True); s[s==0] = 1.0
    Z = X0 / s
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    expl = (S**2) / max(1.0, Z.shape[0]-1); expl = expl / max(1e-12, expl.sum())
    load = Vt.T; importance = (load**2 * expl).sum(axis=1); importance[importance<0]=0
    s2 = importance.sum(); return importance / (s2 if s2>1e-12 else 1.0)

def learn_feature_weights(rows: List[Dict[str,str]], P: SParams, mode: str):
    name_ref, v0 = build_feature_scores_from_row(rows[0], P); mat = [v0]
    for r in rows[1:]:
        _, v = build_feature_scores_from_row(r, P); mat.append(v)
    if np is None:
        # Uniform over non-constant features
        m = len(mat); k = len(name_ref)
        means = [sum(row[j] for row in mat)/m for j in range(k)]
        vars_ = []
        for j in range(k):
            v = sum((row[j]-means[j])**2 for row in mat)/(m-1 if m>1 else 1)
            vars_.append(v)
        keep = [i for i, vv in enumerate(vars_) if vv > 1e-9]
        w = [0.0]*k
        if keep:
            unif = 1.0/len(keep)
            for i in keep: w[i] = unif
        return name_ref, w
    X = np.array(mat, dtype=float)
    var = X.var(axis=0)
    keep_idx = [i for i, vv in enumerate(var) if vv > 1e-9]
    if not keep_idx: return name_ref, [1.0/len(name_ref)]*len(name_ref)
    Xk = X[:, keep_idx]
    wk_core = pca_weights(Xk) if mode == "pca" else entropy_weights(Xk)
    w_full = [0.0]*X.shape[1]
    for j, idx in enumerate(keep_idx): w_full[idx] = float(wk_core[j])
    s2 = sum(w_full); w_full = ([1.0/len(w_full)]*len(w_full)) if s2<=1e-12 else [w/s2 for w in w_full]
    return name_ref, w_full

def load_csv_rows(path: str) -> List[Dict[str,str]]:
    rows: List[Dict[str,str]] = []
    with open(path, newline="") as f:
        rdr = csv.DictReader(f)
        for r in rdr: rows.append(r)
    if not rows: raise RuntimeError("No rows in input CSV.")
    return rows

def write_csv_scored(path: str, rows: List[Dict[str,str]], extra_order: List[str]):
    base_keys = [k for k in rows[0].keys() if not (k.startswith("score_") or k.startswith("weight_") or k.startswith("contrib_") or k in {"score_mode","composite_score","composite_entropy","composite_pca","composite_rank","composite_tiebreak"})]
    fields = list(base_keys) + ["score_mode","composite_score","composite_entropy","composite_pca","composite_rank","composite_tiebreak"] + extra_order
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for r in rows: w.writerow(r)

def rank_borda(scores_matrix: List[List[float]]) -> List[float]:
    """
    Rank-based composite (Borda): for each feature, assign ranks (higher is better),
    convert to [0,1], and average across features.
    """
    import numpy as _np  # local to avoid requiring numpy globally
    S = _np.array(scores_matrix, dtype=float)  # n x d
    n, d = S.shape
    # ranks per column (argsort twice trick for dense ranks)
    ranks = _np.zeros_like(S)
    for j in range(d):
        order = _np.argsort(S[:, j])
        dense = _np.empty_like(order)
        dense[order] = _np.arange(n)  # 0..n-1 (low..high)
        ranks[:, j] = dense
    # normalize to [0,1]
    ranks01 = ranks / max(1.0, (n-1))
    # average across features
    return ranks01.mean(axis=1).tolist()

def deterministic_tiebreak(row: Dict[str,str]) -> float:
    """
    Tiny deterministic epsilon in [0,1), based on sequences/coords, to kill residual ties.
    """
    key = f"{row.get('guide_idx','')}|{row.get('guide_seq','')}|{row.get('fp_start','')}|{row.get('rp_start','')}|{row.get('fp_seq','')}|{row.get('rp_seq','')}"
    h = hashlib.sha1(key.encode("utf-8")).hexdigest()
    # take 8 hex chars -> 32 bits -> 0..1
    v = int(h[:8], 16) / float(0xFFFFFFFF)
    return v

def score_rows(rows: List[Dict[str,str]], P: SParams, mode_entropy: str, mode_pca: str):
    # Build matrix
    feats, _ = build_feature_scores_from_row(rows[0], P)
    score_mat: List[List[float]] = []
    for r in rows:
        names, s = build_feature_scores_from_row(r, P)
        if names != feats:
            raise RuntimeError("Feature ordering mismatch.")
        score_mat.append(s)

    # Entropy composite
    if np is None:
        # fallback uniform weights
        wE = [1.0/len(feats)] * len(feats)
    else:
        X = np.array(score_mat, dtype=float)
        # compute entropy weights on current matrix
        col_sum = X.sum(axis=0)
        col_sum[col_sum == 0] = 1e-12
        Pmat = X / col_sum
        n = max(2, X.shape[0]); k = 1.0 / math.log(n)
        H = -k * (Pmat * np.log(Pmat + 1e-12)).sum(axis=0)
        D = np.maximum(0.0, 1.0 - H)
        sD = D.sum(); wE = (D / sD).tolist() if sD > 1e-12 else [1.0/len(feats)]*len(feats)

    comp_entropy = [sum(si*w for si, w in zip(s, wE)) for s in score_mat]

    # PCA composite
    if np is None:
        wP = [1.0/len(feats)]*len(feats)
    else:
        X = np.array(score_mat, dtype=float)
        X0 = X - X.mean(axis=0, keepdims=True)
        std = X0.std(axis=0, ddof=1, keepdims=True); std[std==0] = 1.0
        Z = X0 / std
        U, S, Vt = np.linalg.svd(Z, full_matrices=False)
        expl = (S**2) / max(1.0, Z.shape[0]-1)
        expl = expl / max(1e-12, expl.sum())
        load = Vt.T
        importance = (load**2 * expl).sum(axis=1)
        importance[importance<0]=0
        sI = importance.sum()
        wP = (importance / sI).tolist() if sI>1e-12 else [1.0/len(feats)]*len(feats)

    comp_pca = [sum(si*w for si, w in zip(s, wP)) for s in score_mat]

    # Rank-based composite
    try:
        comp_rank = rank_borda(score_mat)
    except Exception:
        # safe fallback: normalized sum
        comp_rank = [sum(s)/len(s) for s in score_mat]

    # Normalize each composite to [0,1] to mix fairly
    def norm01(vs: List[float]) -> List[float]:
        vmin = min(vs); vmax = max(vs)
        if vmax - vmin < 1e-12: return [0.5]*len(vs)
        return [(v - vmin)/(vmax - vmin) for v in vs]

    cE = norm01(comp_entropy)
    cP = norm01(comp_pca)
    cR = norm01(comp_rank)

    # Final composite (weights chosen empirically to spread scores)
    final = [0]*len(rows)
    for i in range(len(rows)):
        base = 0.60*cE[i] + 0.35*cP[i] + 0.05*cR[i]
        eps = deterministic_tiebreak(rows[i]) * 1e-6  # micro epsilon for last-resort tie kill
        final[i] = base + eps

    # Build output rows
    out_rows: List[Dict[str,str]] = []
    for i, r in enumerate(rows):
        r_out = dict(r)
        r_out["score_mode"] = f"mix(entropy+pca+rank)"
        r_out["composite_entropy"] = f"{cE[i]:.8f}"
        r_out["composite_pca"]     = f"{cP[i]:.8f}"
        r_out["composite_rank"]    = f"{cR[i]:.8f}"
        r_out["composite_tiebreak"]= f"{deterministic_tiebreak(r):.8f}"
        r_out["composite_score"]   = f"{final[i]:.8f}"
        # also emit per-feature scores & both weight sets so you can inspect
        names, s = build_feature_scores_from_row(r, P)
        for n, si in zip(names, s):
            r_out[f"score_{n}"] = f"{si:.6f}"
        for n, wi in zip(names, wE):
            r_out[f"weightE_{n}"] = f"{wi:.6f}"
        for n, wi in zip(names, wP):
            r_out[f"weightP_{n}"] = f"{wi:.6f}"
        out_rows.append(r_out)

    # Return feature order and out rows
    return feats, out_rows

def scoring_main(args) -> str:
    rows = load_csv_rows(args.input_csv)
    P = SParams(
        gc_primer=parse_range_2floats(args.gc_primer, "gc-primer"),
        amplicon_len=parse_range_2floats(args.amplicon_len, "amplicon-len", lo_ok=1.0, hi_ok=1e9),
        guide_margin_bp=float(args.guide_margin),
        tm_delta_warn=float(args.tm_delta),
        cross_dimer_thresh=float(args.cross_dimer_thresh),
        self_dimer_thresh=float(args.self_dimer_thresh),
        seed_gc_cas9=parse_range_2floats(args.seed_gc_cas9, "seed-gc-cas9"),
        seed_gc_other=parse_range_2floats(args.seed_gc_other, "seed-gc-other"),
    )
    feats, out_rows = score_rows(rows, P, mode_entropy="entropy", mode_pca="pca")

    extra_cols: List[str] = []
    # weights & scores
    for n in feats: extra_cols += [f"weightE_{n}"]
    for n in feats: extra_cols += [f"weightP_{n}"]
    for n in feats: extra_cols += [f"score_{n}"]

    write_csv_scored(args.output_csv, out_rows, extra_cols)

    # Short console summary
    print(f"[OK] Scored {len(out_rows)} rows with mixed composite (entropy+pca+rank+deterministic-ε).")
    # show top 10 by composite_score
    top = sorted(out_rows, key=lambda r: float(r["composite_score"]), reverse=True)[:10]
    print("Top 10 by composite_score:")
    for r in top:
        print(f"  score={r['composite_score']}  guide_idx={r.get('guide_idx')}  fp@{r.get('fp_start')} rp@{r.get('rp_start')}  amp={r.get('amplicon_len')}")

    return args.output_csv

# ============================================================
# CLI
# ============================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate primer/crRNA features with staged prefilter (→ ~1000 diverse) & tie-resistant scoring. Background FASTA optional.")
    sub = p.add_subparsers(dest="cmd", required=True)

    # ---- generate
    g = sub.add_parser("generate", help="Generate features CSV (after heavy stage)")
    g.add_argument("--target-fasta", type=str, required=True)
    g.add_argument("--cas-type", type=validate_cas, required=True)
    g.add_argument("--background-fasta", type=str, default=None, help="Optional background FASTA for BLAST")
    g.add_argument("--top-k", type=int, default=10)
    g.add_argument("--primer-len", type=lambda s: parse_range_2ints(s, "primer-len"), default="30-35")
    g.add_argument("--amplicon-len", type=lambda s: parse_range_2ints(s, "amplicon-len"), default="100-250")
    g.add_argument("--gc-primer", type=lambda s: parse_range_2floats(s, "gc-primer"), default="30-70")
    g.add_argument("--max-homopolymer", type=int, default=5)
    g.add_argument("--crrna-len", type=int, default=None)
    g.add_argument("--gc-crrna", type=lambda s: parse_range_2floats(s, "gc-crrna"), default="35-55")
    g.add_argument("--max-homopolymer-crrna", type=int, default=5)
    g.add_argument("--cross-dimer-thresh", type=float, default=40.0)
    g.add_argument("--self-dimer-thresh", type=float, default=40.0)
    g.add_argument("--tm-delta-warn", type=float, default=10.0)
    g.add_argument("--blast-word-size", type=int, default=7)
    g.add_argument("--blast-evalue", type=float, default=1000.0)
    g.add_argument("--blast-task", type=str, default="blastn-short")
    g.add_argument("--threads", type=int, default=1)
    g.add_argument("--guide-margin", type=int, default=20)
    g.add_argument("--on-target-k", type=int, default=12)
    g.add_argument("--offtarget-max-allow", type=float, default=0.80)
    g.add_argument("--write-csv", type=str, default="pairs_features.csv")
    # staged prefilter controls (default target 1000)
    g.add_argument("--prefilter-target", type=int, default=1000, help="How many pairs to keep before heavy stage")
    g.add_argument("--similarity-bin", type=int, default=5, help="Grid bin (bp) for FP/RP starts to avoid near-duplicates")
    g.add_argument("--per-guide-cap-mult", type=float, default=2.5, help="Max per-guide ≈ base_quota * multiplier")
    g.add_argument("--skip-blast", action="store_true", help="Force skip BLAST even if background provided")
    g.add_argument("--skip-vienna", action="store_true", help="Skip guide secondary structure (fast)")
    g.add_argument("--random-seed", type=int, default=None, help="Optional seed for any randomized steps (not used by default)")

    # ---- score
    s = sub.add_parser("score", help="Score an existing features CSV")
    s.add_argument("--input-csv", required=True)
    s.add_argument("--output-csv", required=True)
    s.add_argument("--gc-primer", default="30-70")
    s.add_argument("--amplicon-len", default="100-250")
    s.add_argument("--guide-margin", type=float, default=20.0)
    s.add_argument("--tm-delta", type=float, default=10.0)
    s.add_argument("--cross-dimer-thresh", type=float, default=40.0)
    s.add_argument("--self-dimer-thresh", type=float, default=40.0)
    s.add_argument("--seed-gc-cas9", default="40-60")
    s.add_argument("--seed-gc-other", default="35-55")

    # ---- gen-score combined
    gs = sub.add_parser("gen-score", help="Generate features and immediately score them")
    # reuse generator args
    for a in g._actions[1:]:
        if a.dest not in {"write_csv"}: gs._add_action(a)
    gs.add_argument("--write-csv", type=str, default="pairs_features2.csv")
    # scoring outputs/params
    gs.add_argument("--scored-csv", type=str, default="pairs_scored.csv")
    gs.add_argument("--seed-gc-cas9", default="40-60")
    gs.add_argument("--seed-gc-other", default="35-55")
    return p

# ============================================================
# Entry
# ============================================================

def main():
    p = build_parser(); args = p.parse_args()
    if args.cmd == "generate":
        out = generator_main(args)
        if not out: sys.exit(2)
    elif args.cmd == "score":
        scoring_main(args)
    elif args.cmd == "gen-score":
        gen_csv = generator_main(args)
        if not gen_csv: sys.exit(2)
        # build scoring namespace
        sargs = argparse.Namespace(
            input_csv=gen_csv,
            output_csv=args.scored_csv,
            gc_primer=f"{args.gc_primer[0]}-{args.gc_primer[1]}",
            amplicon_len=f"{args.amplicon_len[0]}-{args.amplicon_len[1]}",
            guide_margin=float(args.guide_margin),
            tm_delta=float(args.tm_delta_warn),
            cross_dimer_thresh=float(args.cross_dimer_thresh),
            self_dimer_thresh=float(args.self_dimer_thresh),
            seed_gc_cas9=args.seed_gc_cas9,
            seed_gc_other=args.seed_gc_other,
        )
        scoring_main(sargs)
    else:
        p.error("Unknown subcommand.")

if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        print("\nERROR running an external tool (makeblastdb/blastn or RNAfold).", file=sys.stderr)
        print(e, file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr); sys.exit(1)
