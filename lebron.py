
from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Tuple, Optional, List, Dict
import argparse
import os
import sys
import re
import csv
import tempfile
import subprocess
from math import isfinite

def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")
    return seq.translate(table)[::-1]

def gc_pct(seq: str) -> float:
    s = seq.upper()
    g = s.count("G"); c = s.count("C")
    atgc = sum(s.count(x) for x in "ACGT")
    return 100.0 * (g + c) / max(1, atgc)

def has_homopolymer(seq: str, max_run: int) -> bool:
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            run += 1
            if run >= max_run:
                return True
        else:
            run = 1
    return False

def wallace_tm(seq: str) -> float:
    """Very rough Tm (2*(A+T) + 4*(G+C)); OK for relative balance checks."""
    s = seq.upper()
    a = s.count("A"); t = s.count("T"); g = s.count("G"); c = s.count("C")
    return 2.0*(a+t) + 4.0*(g+c)

def gentle_gc_clamp_ok(seq: str) -> int:
    """Return GC count in last 3 nt (want exactly 1 for a gentle clamp)."""
    tail = seq[-3:] if len(seq) >= 3 else seq
    return tail.count("G") + tail.count("C")

def max_complementarity_percent(a: str, b: str) -> float:
    """Max fraction of complementary matches across all shifts (as %)."""
    comp = str.maketrans("ACGT", "TGCA")
    aU = a.upper()
    bC = b.upper().translate(comp)[::-1]
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
        description="Steps 1–7: Config + crRNA + primers + pairing + structure + BLAST specificity"
    )

    # Required inputs
    p.add_argument("--target-fasta", type=str, required=True,
                   help="Path to FASTA with target (first record used)")
    p.add_argument("--cas-type", type=validate_cas, required=True,
                   help="Cas nuclease: cas12a | cas9 | cas13")
    p.add_argument("--background-fasta", type=str, required=True,
                   help="Background genome/transcriptome FASTA for BLAST off-target checks")
    p.add_argument("--top-k", type=int, default=10,
                   help="How many paired sets to print to console (CSV contains all)")

    # Primer & amplicon
    p.add_argument("--primer-len", type=lambda s: parse_range_2ints(s, "primer-len"),
                   default="30-35")
    p.add_argument("--amplicon-len", type=lambda s: parse_range_2ints(s, "amplicon-len"),
                   default="100-250")
    p.add_argument("--gc-primer", type=lambda s: parse_range_2floats(s, "gc-primer"),
                   default="30-70")
    p.add_argument("--max-homopolymer", type=int, default=5)

    # crRNA
    p.add_argument("--crrna-len", type=int, default=None)
    p.add_argument("--gc-crrna", type=lambda s: parse_range_2floats(s, "gc-crrna"),
                   default="35-55")
    p.add_argument("--max-homopolymer-crrna", type=int, default=5)

    # Structure thresholds + ΔTm guidance
    p.add_argument("--cross-dimer-thresh", type=float, default=40.0)
    p.add_argument("--self-dimer-thresh", type=float, default=40.0)
    p.add_argument("--tm-delta-warn", type=float, default=10.0,
                   help="Warn when |Tm(FP)-Tm(RP)| exceeds this (°C)")

    # BLAST
    p.add_argument("--blast-word-size", type=int, default=4)
    p.add_argument("--blast-evalue", type=float, default=1000.0)
    p.add_argument("--blast-task", type=str, default="blastn")
    p.add_argument("--threads", type=int, default=1)

    # Pairing geometry margin
    p.add_argument("--guide-margin", type=int, default=20,
                   help="Guide must be ≥ this many bp from each primer end (default 20)")

    # Output CSV
    p.add_argument("--write-csv", type=str, default="pairs_features.csv",
                   help="Output CSV with raw feature metrics (default pairs_features.csv)")

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
    )

    print("\n[CONFIG]")
    for k, v in asdict(cfg).items():
        print(f"  {k}: {v}")
    print(f"\n[INPUT] target {target_name}: {len(target_seq)} bp")
    print(f"[INPUT] background FASTA: {cfg.background_fasta}\n")

    return cfg, target_name, target_seq, cfg.background_fasta, args.write_csv


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

    # Light filters
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
                # Sanity margins
                left_gap = g.protospacer_start - amplicon_start
                right_gap = amplicon_end - (g.protospacer_start + cfg.crrna_len)
                if left_gap < margin or right_gap < margin:
                    continue
                pairs.append(PairedSet(
                    guide_index=gi, guide=g, fp=fp, rp=rp,
                    amplicon_len=amplicon, metrics={}, flags={}
                ))
    return pairs



def annotate_structure_and_tm(ps: PairedSet, cfg: Config):
    # Tm
    tm_fp = wallace_tm(ps.fp.seq)
    tm_rp = wallace_tm(ps.rp.seq)
    delta_tm = abs(tm_fp - tm_rp)

    # Dimer/hairpin (simple max complementarity)
    fp_self = max_complementarity_percent(ps.fp.seq, ps.fp.seq)
    rp_self = max_complementarity_percent(ps.rp.seq, ps.rp.seq)
    fp_rp_cross = max_complementarity_percent(ps.fp.seq, ps.rp.seq)

    # Flags if exceeding thresholds
    if fp_self > cfg.self_dimer_thresh_pct:
        ps.flags["fp_self_dimer"] = f"{fp_self:.1f}%"
    if rp_self > cfg.self_dimer_thresh_pct:
        ps.flags["rp_self_dimer"] = f"{rp_self:.1f}%"
    if fp_rp_cross > cfg.cross_dimer_thresh_pct:
        ps.flags["fp_rp_cross_dimer"] = f"{fp_rp_cross:.1f}%"
    if delta_tm > cfg.tm_delta_warn:
        ps.flags["delta_tm_high"] = f"{delta_tm:.1f}C"

    # Save raw metrics
    ps.metrics.update({
        "tm_fp_C": tm_fp,
        "tm_rp_C": tm_rp,
        "delta_tm_C": delta_tm,
        "fp_self_dimer_pct": fp_self,
        "rp_self_dimer_pct": rp_self,
        "fp_rp_cross_dimer_pct": fp_rp_cross
    })

# =========================
# Step 6: Guide context features
# =========================

def annotate_context(ps: PairedSet, cfg: Config):
    amp_start = ps.fp.start
    amp_end = ps.rp.start + ps.rp.length  # exclusive
    amp_mid = amp_start + ps.amplicon_len / 2.0

    g_start = ps.guide.protospacer_start
    g_end = g_start + cfg.crrna_len
    g_mid = (g_start + g_end) / 2.0

    left_margin = g_start - amp_start
    right_margin = amp_end - g_end
    centered_abs_bp = abs(g_mid - amp_mid)

    # Seed region (Cas12a 1–8 from 5' end of protospacer; Cas9 seed is PAM-proximal)
    seed_start = 0
    seed_len = 8
    if ps.guide.cas_type == "cas9":
        # PAM-proximal on the protospacer 3' end (positions 12–20 typically);
        # take last 8 bases as a simple proxy
        seed_seq = ps.guide.protospacer_seq[-8:]
    else:
        # Cas12a/Cas13: 5' end seed proxy
        seed_seq = ps.guide.protospacer_seq[seed_start:seed_start+seed_len]
    seed_gc = gc_pct(seed_seq)
    seed_hpoly = has_homopolymer(seed_seq, 4)

    ps.metrics.update({
        "left_margin_bp": float(left_margin),
        "right_margin_bp": float(right_margin),
        "amplicon_center_offset_bp": float(centered_abs_bp),
        "seed_gc_pct": seed_gc,
        "seed_has_homopolymer": 1.0 if seed_hpoly else 0.0
    })

# =========================
# Step 7: BLAST specificity (FP, RP, crRNA)
# =========================

class BlastCache:
    def __init__(self):
        self.cache: Dict[str, Tuple[float, str]] = {}

    def get(self, seq: str):
        return self.cache.get(seq)

    def set(self, seq: str, val: Tuple[float, str]):
        self.cache[seq] = val

def which(cmd: str) -> Optional[str]:
    from shutil import which as _which
    return _which(cmd)

def ensure_blast_db(background_fasta: str, workdir: str) -> str:
    """Create a BLAST DB in workdir if not present. Returns DB prefix path."""
    mk = which("makeblastdb")
    if not mk:
        raise RuntimeError("makeblastdb not found in PATH.")
    db_prefix = os.path.join(workdir, "background_db")
    # Build if needed (we always rebuild in temp workdir)
    cmd = [
        "makeblastdb", "-in", background_fasta, "-dbtype", "nucl",
        "-parse_seqids", "-out", db_prefix
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return db_prefix

def blast_crossreact_score(query_seq: str,
                           db_prefix: str,
                           cfg: Config,
                           workdir: str) -> Tuple[float, str]:
    """
    Run blastn for one query_seq; return (max_score, sseqid_of_max).
    Score = max over hits of (pident/100) * (aln_length / query_length).
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
        "-task", cfg.blast_task,
        "-dust", "no",
        "-soft_masking", "false",
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sstrand",
        "-out", out_tsv
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    max_score = 0.0
    max_id = ""
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    continue
                pident = float(parts[2])
                aln_len = float(parts[3])
                sseqid = parts[1]
                score = (pident/100.0) * (aln_len / max(1.0, float(len(query_seq))))
                if score > max_score:
                    max_score = score
                    max_id = sseqid
    return (max_score, max_id)

def annotate_specificity_with_blast(pairs: List[PairedSet], cfg: Config, background_fasta: str):
    if which("blastn") is None:
        raise RuntimeError("blastn not found in PATH.")
    with tempfile.TemporaryDirectory(prefix="primedrpa_blast_") as td:
        db_prefix = ensure_blast_db(background_fasta, td)
        cache = BlastCache()
        for ps in pairs:
            # FP
            c = cache.get(ps.fp.seq)
            if c is None:
                c = blast_crossreact_score(ps.fp.seq, db_prefix, cfg, td)
                cache.set(ps.fp.seq, c)
            ps.metrics["fp_offtarget_score"] = c[0]
            ps.metrics["fp_offtarget_hit"] = c[1] or ""
            # RP
            c = cache.get(ps.rp.seq)
            if c is None:
                c = blast_crossreact_score(ps.rp.seq, db_prefix, cfg, td)
                cache.set(ps.rp.seq, c)
            ps.metrics["rp_offtarget_score"] = c[0]
            ps.metrics["rp_offtarget_hit"] = c[1] or ""
            # crRNA
            gseq = ps.guide.protospacer_seq
            c = cache.get(gseq)
            if c is None:
                c = blast_crossreact_score(gseq, db_prefix, cfg, td)
                cache.set(gseq, c)
            ps.metrics["crrna_offtarget_score"] = c[0]
            ps.metrics["crrna_offtarget_hit"] = c[1] or ""

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
        "fp_start", "fp_len", "fp_gc_pct", "fp_seq",
        "rp_start", "rp_len", "rp_gc_pct", "rp_seq",
        "amplicon_len",
        "left_margin_bp", "right_margin_bp", "amplicon_center_offset_bp",
        # thermo/structure
        "tm_fp_C", "tm_rp_C", "delta_tm_C",
        "fp_self_dimer_pct", "rp_self_dimer_pct", "fp_rp_cross_dimer_pct",
        # seed/context
        "seed_gc_pct", "seed_has_homopolymer",
        # specificity
        "fp_offtarget_score", "fp_offtarget_hit",
        "rp_offtarget_score", "rp_offtarget_hit",
        "crrna_offtarget_score", "crrna_offtarget_hit",
        # flags (semicolon-joined)
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
                "rp_start": ps.rp.start,
                "rp_len": ps.rp.length,
                "rp_gc_pct": f"{ps.rp.gc_pct:.1f}",
                "rp_seq": ps.rp.seq,
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
                "seed_gc_pct": f"{ps.metrics.get('seed_gc_pct', 0):.1f}",
                "seed_has_homopolymer": int(ps.metrics.get('seed_has_homopolymer', 0)),
                "fp_offtarget_score": f"{ps.metrics.get('fp_offtarget_score', 0):.3f}",
                "fp_offtarget_hit": ps.metrics.get('fp_offtarget_hit', ""),
                "rp_offtarget_score": f"{ps.metrics.get('rp_offtarget_score', 0):.3f}",
                "rp_offtarget_hit": ps.metrics.get('rp_offtarget_hit', ""),
                "crrna_offtarget_score": f"{ps.metrics.get('crrna_offtarget_score', 0):.3f}",
                "crrna_offtarget_hit": ps.metrics.get('crrna_offtarget_hit', ""),
                "flags": ";".join([f"{k}={v}" for k, v in ps.flags.items()]) if ps.flags else ""
            }
            w.writerow(row)
    print(f"[WRITE] CSV → {out_path}  ({len(pairs)} rows)")

# =========================
# Main (Steps 1–7)
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

        # Step 5 + Step 6 per pair
        for ps in pairs:
            annotate_structure_and_tm(ps, cfg)
            annotate_context(ps, cfg)

        # Step 7 (BLAST specificity)
        print("[BLAST] Building DB and checking off-targets for FP, RP, and crRNA...")
        annotate_specificity_with_blast(pairs, cfg, bg)

        # Console preview of top_k
        show_n = min(cfg.top_k, len(pairs))
        for ps in pairs[:show_n]:
            g = ps.guide
            print(f"\n[SET] GuideIdx={ps.guide_index} {g.cas_type}{g.strand} "
                  f"{'PAM='+g.pam_seq if g.pam_seq else 'PAM=NA'} prot@{g.protospacer_start} len={len(g.protospacer_seq)} GC={g.gc_pct:.1f}%")
            print(f"  FP  start={ps.fp.start} len={ps.fp.length} GC={ps.fp.gc_pct:.1f}%  Tm={ps.metrics['tm_fp_C']:.1f}°C  self%={ps.metrics['fp_self_dimer_pct']:.1f}  seq={ps.fp.seq}")
            print(f"  RP  start={ps.rp.start} len={ps.rp.length} GC={ps.rp.gc_pct:.1f}%  Tm={ps.metrics['tm_rp_C']:.1f}°C  self%={ps.metrics['rp_self_dimer_pct']:.1f}  seq={ps.rp.seq}")
            print(f"  Cross-dimer%={ps.metrics['fp_rp_cross_dimer_pct']:.1f}   ΔTm={ps.metrics['delta_tm_C']:.1f}°C")
            lm = ps.metrics['left_margin_bp']; rm = ps.metrics['right_margin_bp']
            print(f"  Amplicon={ps.amplicon_len} bp   margins L={lm:.0f} R={rm:.0f}   center_offset={ps.metrics['amplicon_center_offset_bp']:.1f} bp")
            print(f"  Seed GC={ps.metrics['seed_gc_pct']:.1f}%   seed_homopolymer={bool(ps.metrics['seed_has_homopolymer'])}")
            print(f"  Off-target FP={ps.metrics['fp_offtarget_score']:.3f}  RP={ps.metrics['rp_offtarget_score']:.3f}  crRNA={ps.metrics['crrna_offtarget_score']:.3f}")
            if ps.flags:
                print(f"  Flags: {ps.flags}")

        # CSV
        write_pairs_csv(pairs, out_csv, tname)
        print("\n[OK] Steps 1–7 complete. CSV contains raw per-feature metrics for each pair.")
    except subprocess.CalledProcessError as e:
        print("\nERROR running an external tool (makeblastdb/blastn).", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)
