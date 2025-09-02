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
    background_fasta: str

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

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Config + crRNA + primers + pairing + features + BLAST with mismatches (per-component only)"
    )
    # Required
    p.add_argument("--target-fasta", type=str, required=True)
    p.add_argument("--cas-type", type=validate_cas, required=True)
    p.add_argument("--background-fasta", type=str, required=True)
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
    return p

def step1_get_config_and_sequences(argv: Optional[List[str]] = None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    target_name, target_seq = load_first_fasta_seq(args.target_fasta)

    if not os.path.exists(args.background_fasta):
        parser.error(f"--background-fasta not found: {args.background_fasta}")

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
        background_fasta=os.path.abspath(args.background_fasta),

        guide_margin_bp=max(0, int(args.guide_margin)),
        on_target_k=max(6, int(args.on_target_k)),
        offtarget_max_allow=float(args.offtarget_max_allow),
    )

    print("\n[CONFIG]")
    for k, v in asdict(cfg).items():
        print(f"  {k}: {v}")
    print(f"\n[INPUT] target {target_name}: {len(target_seq)} bp")
    print(f"[INPUT] background FASTA: {cfg.background_fasta}\n")

    return cfg, target_name, target_seq, cfg.background_fasta, args.write_csv

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
        # FP must end before guide start - margin
        left_ok = [fp for fp in fps if (fp.start + fp.length) <= (g.protospacer_start - margin)]
        # RP must start after guide end + margin
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
# Feature annotations (thermo, context, cross-talk, on-target)
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
    t_rc = revcomp(t)
    k = min(cfg.on_target_k, ps.fp.length, ps.rp.length)

    rp_site_top = revcomp(ps.rp.seq)

    fp_hits = len(find_all_indices(t, ps.fp.seq)) + len(find_all_indices(t, revcomp(ps.fp.seq)))
    rp_hits = len(find_all_indices(t, rp_site_top)) + len(find_all_indices(t, revcomp(rp_site_top)))

    fp_k = ps.fp.seq[-k:]
    rp_k = rp_site_top[-k:]
    fp_k_hits = len(find_all_indices(t, fp_k))
    rp_k_hits = len(find_all_indices(t, rp_k))

    fp_mm1 = count_approx_matches(t, ps.fp.seq, 1, extra_string=revcomp(ps.fp.seq))
    rp_mm1 = count_approx_matches(t, rp_site_top, 1, extra_string=revcomp(rp_site_top))
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
        "fp_full_hits_mm1_both": float(fp_mm1),
        "rp_full_hits_mm1_both": float(rp_mm1),

        "fp_3p_k": float(k),
        "fp_3p_k_hits_target": float(fp_k_hits),
        "rp_3p_k_hits_target": float(rp_k_hits),
        "fp_3p_k_hits_mm1": float(fp_k_mm1),
        "rp_3p_k_hits_mm1": float(rp_k_mm1),

        "amplicon_hits_on_target": float(amplicon_hits),
    })

    if fp_k_hits > 1:
        ps.flags["fp_3p_k_nonunique"] = f"{fp_k_hits}"
    if rp_k_hits > 1:
        ps.flags["rp_3p_k_nonunique"] = f"{rp_k_hits}"
    if amplicon_hits != 1:
        ps.flags["amplicon_not_unique"] = str(amplicon_hits)

# =========================
# ViennaRNA / RNAfold integration
# =========================

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
# BLAST (mismatch-aware) — PER-COMPONENT ONLY
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
    """
    Return ALL hits as dicts:
      {'sseqid','pident','aln_len','mismatch','qstart','qend','sstart','send','sstrand','score','is_exact_full'}
    score = (pident/100) * (aln_len / qlen)
    """
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
        "-task", cfg.blast_task,       # consider 'blastn-short' for short oligos
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
        return None, None, 0, ""  # top_any, top_imp, n_exact, example_exact_id
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
                # exact full-length presence
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

# =========================
# CSV writer
# =========================

def write_pairs_csv(pairs: List[PairedSet], out_path: str, target_name: str):
    fields = [
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
        # on-target uniqueness
        "fp_full_hits_target_both", "rp_full_hits_target_both",
        "fp_full_hits_mm1_both", "rp_full_hits_mm1_both",
        "fp_3p_k", "fp_3p_k_hits_target", "rp_3p_k_hits_target",
        "fp_3p_k_hits_mm1", "rp_3p_k_hits_mm1",
        "amplicon_hits_on_target",
        # mismatch-aware individual off-targets (best-any + best-imperfect) + exact flags
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
        # flags
        "flags"
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for ps in pairs:
            g = ps.guide
            row = {
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
                "fp_full_hits_mm1_both": int(ps.metrics.get('fp_full_hits_mm1_both', 0)),
                "rp_full_hits_mm1_both": int(ps.metrics.get('rp_full_hits_mm1_both', 0)),
                "fp_3p_k": int(ps.metrics.get("fp_3p_k", 0)),
                "fp_3p_k_hits_target": int(ps.metrics.get('fp_3p_k_hits_target', 0)),
                "rp_3p_k_hits_target": int(ps.metrics.get('rp_3p_k_hits_target', 0)),
                "fp_3p_k_hits_mm1": int(ps.metrics.get('fp_3p_k_hits_mm1', 0)),
                "rp_3p_k_hits_mm1": int(ps.metrics.get('rp_3p_k_hits_mm1', 0)),
                "amplicon_hits_on_target": int(ps.metrics.get('amplicon_hits_on_target', 0)),

                # per-component off-targets
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

# =========================
# Main
# =========================

if __name__ == "__main__":
    try:
        cfg, tname, tseq, bg, out_csv = step1_get_config_and_sequences()

        # Step 2
        guides = step2_discover_crrnas(tseq, cfg)
        if not guides:
            print("\n[NOTE] No crRNA candidates passed light filtering.")
            sys.exit(0)

        # Step 3
        fps = candidate_primers_plus_strand(tseq, cfg)
        rps = candidate_primers_minus_strand(tseq, cfg)
        print(f"\n[PRIMERS] {len(fps)} forward; {len(rps)} reverse candidates passed cheap filters "
              f"(GC {cfg.gc_primer_range[0]}–{cfg.gc_primer_range[1]}%, homopolymer<{cfg.max_homopolymer}, 3' clamp=1 GC).")
        if not fps or not rps:
            print("[NOTE] No primer candidates after cheap filters.")
            sys.exit(0)

        # Step 4
        pairs = pair_primers_around_guides(tseq, guides, fps, rps, cfg)
        print(f"[PAIRING] {len(pairs)} FP/RP pairs span guides with margin {cfg.guide_margin_bp} bp and amplicon "
              f"{cfg.amplicon_len_range[0]}–{cfg.amplicon_len_range[1]} bp.")
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

        # Step 7: BLAST specificity (per-component only)
        print("[BLAST] Building DB and checking off-targets for FP, RP, and crRNA (per-component, mismatch-aware)...")
        annotate_specificity_with_blast(pairs, cfg, bg)

        # Console preview
        show_n = min(cfg.top_k, len(pairs))
        for ps in pairs[:show_n]:
            g = ps.guide
            print(f"\n[SET] GuideIdx={ps.guide_index} {g.cas_type}{g.strand} "
                  f"{'PAM='+g.pam_seq if g.pam_seq else 'PAM=NA'} prot@{g.protospacer_start} len={len(g.protospacer_seq)} GC={g.gc_pct:.1f}%")
            print(f"  FP  start={ps.fp.start} len={ps.fp.length} GC={ps.fp.gc_pct:.1f}%  Tm={ps.metrics['tm_fp_C']:.1f}°C  self%={ps.metrics['fp_self_dimer_pct']:.1f}  seq={ps.fp.seq}")
            print(f"  RP  start={ps.rp.start} len={ps.rp.length} GC={ps.rp.gc_pct:.1f}%  Tm={ps.metrics['tm_rp_C']:.1f}°C  self%={ps.metrics['rp_self_dimer_pct']:.1f}  seq={ps.rp.seq}")
            print(f"  Cross-dimer%={ps.metrics['fp_rp_cross_dimer_pct']:.1f}   ΔTm={ps.metrics['delta_tm_C']:.1f}°C   3'cross_run={int(ps.metrics['fp_rp_3p_cross_run'])}")
            lm = ps.metrics['left_margin_bp']; rm = ps.metrics['right_margin_bp']
            print(f"  Amplicon={ps.amplicon_len} bp   margins L={lm:.0f} R={rm:.0f}   center_offset={ps.metrics['amplicon_center_offset_bp']:.1f} bp  GC={ps.metrics['amplicon_gc_pct']:.1f}%")
            print(f"  Seed seq={ps.metrics['seed_seq']}  GC={ps.metrics['seed_gc_pct']:.1f}%  max_run={int(ps.metrics['seed_max_run'])}")

            # Per-component off-target summary (exact match presence retained)
            def _sumline(prefix: str, label: str):
                exact = "YES" if ps.metrics.get(f"{prefix}_has_exact_full_hit", 0) else "no"
                mm_any = ps.metrics.get(f"{prefix}_offtarget_mm", 0)
                mm_imp = ps.metrics.get(f"{prefix}_offtarget_mm_imp", 0)
                src = ps.metrics.get(f"{prefix}_offtarget_from","")
                hit = ps.metrics.get(f"{prefix}_offtarget_hit","")
                print(f"  {label}: top={hit or '-'}  mm_any={int(mm_any)}  mm_imp={int(mm_imp)} ({src})  exact_full={exact}")

            _sumline("fp", "FP off-target")
            _sumline("rp", "RP off-target")
            _sumline("crrna", "crRNA off-target")

            if ps.flags:
                print(f"  Flags: {ps.flags}")

        # CSV
        write_pairs_csv(pairs, out_csv, tname)
        print("\n[OK] Features + mismatch-aware per-component specificity complete. CSV contains exact-match flags and mismatch counts.")
    except subprocess.CalledProcessError as e:
        print("\nERROR running an external tool (makeblastdb/blastn or RNAfold).", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)
