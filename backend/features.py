import pandas as pd
import primer3
import ViennaRNA
import RNA
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
class RankingFilter:
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

    CATEGORY_FEATURES = {
        "Off-targets": ["crrna_offtarget_mm_imp", "fp_offtarget_mm_imp", "rp_offtarget_mm_imp"],
        "3′ end stability": ["fp_3p_self_run", "rp_3p_self_run", "fp_rp_3p_cross_run"],
        "Primer dimers / secondary structure": ["fp_self_dimer_pct", "rp_self_dimer_pct", "fp_rp_cross_dimer_pct"],
        "Guide seed & structure": ["seed_len", "seed_max_run", "seed_gc_pct", "guide_mfe_kcal", "guide_seed_unpaired_frac"],
        "Tm balance": ["delta_tm_C"],
        "Cross-interactions": ["overlap_protospacer", "overlap_pam"],
        "Conservation": ["guide_conservation", "fp_conservation", "rp_conservation"],
    }

    def __init__(self, csv_filename):
        self.csv_filename = csv_filename
        self.df = pd.read_csv(csv_filename)

    def filter_and_rank(self):
        df = self.df[
            (self.df["overlap_pam"] == 0) &
            (self.df["overlap_protospacer"] == 0) &
            (self.df["amplicon_len"].between(100, 220)) &
            (self.df["fp_len"].between(28, 36)) &
            (self.df["rp_len"].between(28, 36)) &
            (self.df["fp_gc_pct"].between(35, 65)) &
            (self.df["rp_gc_pct"].between(35, 65)) &
            (self.df["amplicon_gc_pct"].between(35, 65)) &
            (self.df["delta_tm_C"] <= 3.0) &
            (self.df["fp_3p_self_run"] <= 2) &
            (self.df["rp_3p_self_run"] <= 2) &
            (self.df["fp_gquad_3p"] == 0) &
            (self.df["rp_gquad_3p"] == 0) &
            (self.df["fp_max_hpoly_run"] <= 4) &
            (self.df["rp_max_hpoly_run"] <= 4) &
            (self.df["amplicon_max_run"] <= 5) &
            (self.df["fp_tail_gc"].between(1, 3)) &
            (self.df["rp_tail_gc"].between(1, 3))
        ].copy()

        # Normalize features
        for feature in self.FEATURE_WEIGHTS:
            if feature in df.columns:
                min_val = df[feature].min()
                max_val = df[feature].max()
                if max_val > min_val:
                    df[feature + "_norm"] = (df[feature] - min_val) / (max_val - min_val)
                else:
                    df[feature + "_norm"] = 0

        # Calculate ranking score
        df["ranking_score"] = 0
        for feature, weight in self.FEATURE_WEIGHTS.items():
            norm_feature = feature + "_norm"
            if norm_feature in df.columns:
                df["ranking_score"] += df[norm_feature] * weight

        max_score = df["ranking_score"].max()
        if max_score > 0:
            df["ranking_score_normalized"] = df["ranking_score"] / max_score
        else:
            df["ranking_score_normalized"] = df["ranking_score"]
        df = df.sort_values(by="ranking_score_normalized", ascending=False)
        df.to_csv("pairs_ranked.csv", index=False)

        # Calculate subscores for each category
        for category, features in self.CATEGORY_FEATURES.items():
            total_weight = sum(self.FEATURE_WEIGHTS[f] for f in features)
            def calc_subscore(row):
                return sum(row[f + "_norm"] * self.FEATURE_WEIGHTS[f] for f in features if f + "_norm" in row) / total_weight
            df[category + "_subscore"] = df.apply(calc_subscore, axis=1)

        print(df[[cat + "_subscore" for cat in self.CATEGORY_FEATURES.keys()]].head())
        return df
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
        rna = _to_rna(guide_seq)
        _, mfe = RNA.fold(rna)  # returns (dotbracket, mfe)
        return float(mfe)

    def guide_seed_unpaired_frac(guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> float:
        """
        Average probability that each base in the seed is unpaired (0..1),
        using ViennaRNA partition function. Higher = more accessible.
        """
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
    def calc_overlap_protospacer(self, guide_start, guide_len, fp_start, fp_len, rp_start, rp_len):
        """
        Returns 1 if protospacer does NOT overlap with either primer (good), else 0 (bad).
        """
        guide_end = guide_start + guide_len - 1
        fp_end = fp_start + fp_len - 1
        rp_end = rp_start + rp_len - 1

        # Check overlap with forward primer
        overlap_fp = not (guide_end < fp_start or guide_start > fp_end)
        # Check overlap with reverse primer
        overlap_rp = not (guide_end < rp_start or guide_start > rp_end)

        return int(not (overlap_fp or overlap_rp))

    def calc_overlap_pam(self, guide_pam_start, pam_len, fp_start, fp_len, rp_start, rp_len):
        """
        Returns 1 if PAM does NOT overlap with either primer (good), else 0 (bad).
        """
        pam_end = guide_pam_start + pam_len - 1
        fp_end = fp_start + fp_len - 1
        rp_end = rp_start + rp_len - 1

        # Check overlap with forward primer
        overlap_fp = not (pam_end < fp_start or guide_pam_start > fp_end)
        # Check overlap with reverse primer
        overlap_rp = not (pam_end < rp_start or guide_pam_start > rp_end)

        return int(not (overlap_fp or overlap_rp))
        


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