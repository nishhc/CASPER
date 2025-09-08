# ...existing code...

import pandas as pd
import primer3
import ViennaRNA
ALPHA_CRISPR = 4.0    # off-target exponent (crRNA)
P_CRISPR     = 1.2
ALPHA_PRIMER = 3.0    # off-target exponent (FP/RP)
P_PRIMER     = 1.0

# Thermo assumptions (tune to your kit/buffer)
DNA_Na_mM = 50.0
DNA_Mg_mM = 8.0
DNA_dNTP_mM = 0.6
DNA_oligo_nM = 250.0

# Seed window for Cas12a (PAM-proximal)
CAS12A_SEED_LEN = 8
class FeatureCalculator:
    def __init__(self, csv_filename):
        self.csv_filename = csv_filename
        self.df = pd.read_csv(csv_filename)
    
    def calc_tm(self, seq: str) -> float:
        return primer3.calcTm(
            seq.upper(),
            mv_conc=DNA_Na_mM, dv_conc=DNA_Mg_mM,
            dntp_conc=DNA_dNTP_mM,
            dna_conc=DNA_oligo_nM
        )
    def percent_mismatch(seq1, seq2):
        # Compare two sequences, return percent mismatches
        min_len = min(len(seq1), len(seq2))
        mismatches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a != b)
        return 100 * mismatches / min_len if min_len > 0 else 0

    def calc_fp_self_dimer_pct(self, fp):
        # Compare forward primer to itself (reverse complement for true dimerization)
        return self.percent_mismatch(fp, fp[::-1])

    def calc_rp_self_dimer_pct(self, rp):
        # Compare reverse primer to itself (reverse complement for true dimerization)
        return self.percent_mismatch(rp, rp[::-1])

    def calc_fp_rp_cross_dimer_pct(self, fp, rp):
        # Compare forward primer to reverse primer
        return self.percent_mismatch(fp, rp)

    def calc_fp_crrna_dimer_pct(self, fp, crrna):
        # Compare forward primer to crRNA
        return self.percent_mismatch(fp, crrna)

    def calc_rp_crrna_dimer_pct(self, rp, crrna):
        # Compare reverse primer to crRNA
        return self.percent_mismatch(rp, crrna)
    def three_prime_self_run(seq: str) -> int:
        """
        Length of the homopolymer run at the 3' end.
        e.g., ...AAAT -> 1 (T), ...AAATT -> 2 (TT), ...AAA -> 3 (AAA)
        """
        s = seq.upper().replace("U","T")
        if not s: return 0
        last = s[-1]
        run = 0
        for ch in reversed(s):
            if ch == last:
                run += 1
            else:
                break
        return run

    def revcomp(seq: str) -> str:
        t = str.maketrans("ACGTUacgtu", "TGCAA tgcaa".replace(" ", ""))
        return seq.translate(t)[::-1]

    def three_prime_cross_run(fp: str, rp: str, tail_len: int = 8) -> int:
        """
        Max contiguous complementarity between the last `tail_len` nt of FP (3' tail)
        and the last `tail_len` nt of RP, after reverse-complementing RP.
        Returns an integer run length (0..tail_len).
        """
        fp_tail = fp[-tail_len:].upper().replace("U","T")
        rp_tail_rc = revcomp(rp[-tail_len:])
        # longest suffix of fp_tail that equals prefix of rp_tail_rc
        best = 0
        for k in range(1, tail_len + 1):
            if fp_tail[-k:] == rp_tail_rc[:k]:
                best = k
        return best
    

    def _to_rna(seq: str) -> str:
        return seq.upper().replace("T", "U")

    def gc_pct(seq: str) -> float:
        s = seq.upper().replace("U","T")
        return (100.0 * (s.count("G") + s.count("C")) / len(s)) if s else 0.0

    def seed_len_max(guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> int:
        """
        Length of the PAM-proximal contiguous seed we consider (cap at guide length).
        For Cas12a, default 8 nt.
        """
        return min(seed_len, len(guide_seq))

    def seed_max_run(guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> int:
        """
        Longest homopolymer run length within the seed window (first seed_len nt).
        """
        s = guide_seq[:seed_len].upper()
        if not s:
            return 0
        runs = [len(r) for r in re.findall(r'(A+|C+|G+|T+|U+)', s)]
        return max(runs) if runs else 0

    def seed_gc_pct(guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> float:
        """
        GC percentage within the seed window (first seed_len nt).
        """
        s = guide_seq[:seed_len]
        return gc_pct(s)

    def guide_mfe_kcal(guide_seq: str) -> float:
        """
        Minimum free energy (kcal/mol) of the guide spacer as RNA.
        More negative = more structured (potentially less accessible).
        """
        if not _HAS_VIENNA:
            # Fallback if ViennaRNA isn't installed
            return float("nan")
        rna = _to_rna(guide_seq)
        _, mfe = RNA.fold(rna)  # returns (dotbracket, mfe)
        return float(mfe)

    def guide_seed_unpaired_frac(guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> float:
        """
        Average probability that each base in the seed is unpaired (0..1),
        using ViennaRNA partition function. Higher = more accessible.
        """
        if not _HAS_VIENNA:
            # Reasonable neutral fallback if ViennaRNA isn't present
            return 0.5
        rna = _to_rna(guide_seq)
        fc = RNA.fold_compound(rna)
        fc.pf()  # compute partition function
        n = len(rna)
        L = min(seed_len, n)
        if L == 0:
            return 0.0
        unpaired_probs = []
        for i in range(1, L + 1):  # 1-based indexing in ViennaRNA
            paired_prob = 0.0
            for j in range(1, n + 1):
                if j == i:
                    continue
                paired_prob += fc.bpp(i, j)
            unpaired_probs.append(max(0.0, 1.0 - paired_prob))
        # clamp to [0,1] and average
        vals = [min(1.0, max(0.0, p)) for p in unpaired_probs]
        return float(sum(vals) / len(vals))
        


FEATURE_WEIGHTS = {

    "crrna_offtarget_mm_imp": 0.12,
    "fp_offtarget_mm_imp": 0.11,
    "rp_offtarget_mm_imp": 0.11,


    "fp_3p_self_run": 0.06,
    "rp_3p_self_run": 0.06,
    "fp_rp_3p_cross_run": 0.06,


    "fp_self_dimer_pct": 0.05,
    "rp_self_dimer_pct": 0.05,
    "fp_rp_cross_dimer_pct": 0.05,


    "seed_len": 0.06,
    "seed_max_run": 0.04,
    "seed_gc_pct": 0.03,
    "guide_mfe_kcal": 0.02,
    "guide_seed_unpaired_frac": 0.02,


    "delta_tm_C": 0.08,


    "overlap_protospacer": 0.02,
    "overlap_pam": 0.02,


    "guide_conservation": 0.02,
    "fp_conservation": 0.005,
    "rp_conservation": 0.005,
}