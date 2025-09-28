# features_only.py
import re
import math
import pandas as pd


try:
    import primer3
except Exception:
    primer3 = None

try:
    import RNA  # ViennaRNA Python bindings
except Exception:
    RNA = None


# ==== Chemistry params (tune as needed) ====
DNA_Na_mM = 50.0
DNA_Mg_mM = 8.0
DNA_dNTP_mM = 0.6
DNA_oligo_nM = 250.0

CAS12A_SEED_LEN = 8


class FeatureCalculator:
    """
    Compute CRISPR/RPA feature columns and append them to the input dataframe.

    Expected (optional) columns in CSV:
      - fp, rp, crrna, amplicon
      - background_seqs (semicolon/pipe/comma-separated DNA strings)
      - guide_aln, fp_aln, rp_aln (semicolon/pipe/comma-separated aligned strings)
      - guide_start, guide_len, fp_start, fp_len, rp_start, rp_len, pam_start, pam_len (ints, 0-based)

    Missing fields are handled; unavailable features default to sensible values/NaNs.
    """

    # ----------------- helpers -----------------
    @staticmethod
    def _revcomp(seq: str) -> str:
        if not isinstance(seq, str):
            return ""
        tbl = str.maketrans("ACGTUacgtu", "TGCAA" + "tgcaa")
        return seq.translate(tbl)[::-1]

    @staticmethod
    def _to_rna(seq: str) -> str:
        return (seq or "").upper().replace("T", "U")

    @staticmethod
    def _gc_pct(seq: str) -> float:
        s = (seq or "").upper().replace("U", "T")
        return (100.0 * (s.count("G") + s.count("C")) / len(s)) if s else 0.0
    @staticmethod
    def _percent_mismatch(seq1: str, seq2: str) -> float:
        s1 = (seq1 or "")
        s2 = (seq2 or "")
        n = min(len(s1), len(s2))
        if n == 0:
            return 0.0
        mm = sum(a != b for a, b in zip(s1[:n], s2[:n]))
        return 100.0 * mm / n

    @staticmethod
    def _three_prime_self_run(seq: str) -> int:
        s = seq.replace("U", "T")
        if not s:
            return 0
        last = s[-1]
        run = 0
        for ch in reversed(s):
            if ch == last:
                run += 1
            else:
                break
        return run

    def _three_prime_cross_run(self, fp: str, rp: str, tail_len: int = 8) -> int:
        fp_tail = (fp or "")[-tail_len:].upper().replace("U", "T")
        rp_tail_rc = self._revcomp((rp or "")[-tail_len:])
        best = 0
        for k in range(1, tail_len + 1):
            if fp_tail[-k:] == rp_tail_rc[:k]:
                best = k
        return best

    @staticmethod
    def _calc_tm(seq: str) -> float:
        s = (seq or "").upper()
        if not s:
            return float("nan")
        try:
            return float(
                primer3.calc_tm(
                    s,
                    mv_conc=DNA_Na_mM,
                    dv_conc=DNA_Mg_mM,
                    dntp_conc=DNA_dNTP_mM,
                    dna_conc=DNA_oligo_nM,
                )
            )
        except Exception:
            return float("nan")

    @staticmethod
    def _seed_len_max(guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> int:
        return min(seed_len, len(guide_seq or ""))

    @staticmethod
    def _seed_max_run(guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> int:
        s = (guide_seq or "")[:seed_len].upper()
        if not s:
            return 0
        runs = [len(r) for r in re.findall(r"(A+|C+|G+|T+|U+)", s)]
        return max(runs) if runs else 0

    def _guide_mfe_kcal(self, guide_seq: str) -> float:
        if not guide_seq or RNA is None:
            return float("nan")
        try:
            rna = self._to_rna(guide_seq)
            _, mfe = RNA.fold(rna)
            return float(mfe)
        except Exception:
            return float("nan")

    def _seed_unpaired_frac(self, guide_seq: str, seed_len: int = CAS12A_SEED_LEN) -> float:
        if not guide_seq or RNA is None:
            return float("nan")
        try:
            rna = self._to_rna(guide_seq)
            fc = RNA.fold_compound(rna)
            fc.pf()
            n = len(rna)
            L = min(seed_len, n)
            if L == 0:
                return float("nan")
            unpaired = []
            # ViennaRNA is 1-based
            for i in range(1, L + 1):
                paired_prob = sum(fc.bpp(i, j) for j in range(1, n + 1) if j != i)
                unpaired.append(max(0.0, 1.0 - paired_prob))
            vals = [min(1.0, max(0.0, p)) for p in unpaired]
            return float(sum(vals) / len(vals)) if vals else float("nan")
        except Exception:
            return float("nan")

    @staticmethod
    def _calc_overlap_interval(a_start, a_len, b_start, b_len) -> bool:
        """Return True if [a] overlaps [b]; tolerate None."""
        try:
            a0, a1 = int(a_start), int(a_start) + int(a_len) - 1
            b0, b1 = int(b_start), int(b_start) + int(b_len) - 1
        except Exception:
            return False
        return not (a1 < b0 or a0 > b1)

    @staticmethod
    def _safe_int(x, default=None):
        try:
            return int(x)
        except Exception:
            return default

    @staticmethod
    def _split_list(cell: str):
        """Split semicolon/pipe/comma-separated cell into list, stripping spaces."""
        if cell is None or (isinstance(cell, float) and math.isnan(cell)):
            return []
        s = str(cell)
        for sep in [";", "|", ","]:
            if sep in s:
                return [t.strip() for t in s.split(sep) if t.strip()]
        s = s.strip()
        return [s] if s else []

    # ----------------- ctor -----------------
    def __init__(self, csv_filename: str):
        self.csv_filename = csv_filename
        self.df = pd.read_csv(csv_filename)

    # ----------------- public API -----------------
    def compute_all_features(self):
        df = self.df

        # Ensure columns exist
        for col in [
            "background_seqs",
        ]:
            if col not in df.columns:
                df[col] = "ATACCAGCTTATTCAATTGGTTGGTCTGGTTGGCCCGTGTGTCATTACGGGTTGGATAAGATAGTAAGTGCAATCTGGGGGTTGTGTTTGAGCGGCGTTTCAGTTGTTTATTTCCCTTTGTTATTCCCTTTGGGGTTGTTGTTTGGTTGTGTGTTTATACCAGCTTATTCAATTCACTTGGTGGTGGTGGCGGGATGGGATGGGTTGGGTTTGTAGATAGTAAGTGCAATCT"

        # Parse list-like cells once
        bg_lists = df["background_seqs"].map(self._split_list)

        # Tm features
        df["fp_tm_C"] = df["forward_primer"].map(self._calc_tm)
        df["rp_tm_C"] = df["backward_primer"].map(self._calc_tm)
        df["delta_tm_C"] = (df["fp_tm_C"] - df["rp_tm_C"]).abs()

        # Dimer proxies (percent mismatches)
        df["fp_self_dimer_pct"] = df["forward_primer"].map(lambda s: self._percent_mismatch(s, self._revcomp(s)))
        df["rp_self_dimer_pct"] = df["backward_primer"].map(lambda s: self._percent_mismatch(s, self._revcomp(s)))
        df["fp_rp_cross_dimer_pct"] = df.apply(
            lambda r: FeatureCalculator._percent_mismatch(r.get("forward_primer"), r.get("backward_primer")), axis=1
        )


        # 3' run features
        # 3' run features
        df["fp_3p_self_run"] = df["forward_primer"].map(lambda s: self._three_prime_self_run(s or ""))
        df["rp_3p_self_run"] = df["backward_primer"].map(lambda s: self._three_prime_self_run(s or ""))

        df["fp_rp_3p_cross_run"] = df.apply(
            lambda r: self._three_prime_cross_run(r.get("forward_primer", ""), r.get("backward_primer", ""), 8), axis=1
        )

        # Seed features (use crRNA as the guide)
        df["seed_len"] = df["crrna"].map(lambda s: self._seed_len_max(s, CAS12A_SEED_LEN))
        df["seed_max_run"] = df["crrna"].map(lambda s: self._seed_max_run(s, CAS12A_SEED_LEN))
        df["seed_gc_pct"] = df["crrna"].map(lambda s: self._gc_pct((s or "")[:CAS12A_SEED_LEN]))

        # Guide structure/accessibility
        df["guide_mfe_kcal"] = df["crrna"].map(self._guide_mfe_kcal)

        # Overlap features (1 = good/no overlap, 0 = bad/overlap)
        def _overlap_protospacer_row(r):
            g_start = self._safe_int(r.get("guide_start"))
            g_len = self._safe_int(r.get("guide_len"))
            fp_s = self._safe_int(r.get("fp_start"))
            fp_l = self._safe_int(r.get("fp_len"))
            rp_s = self._safe_int(r.get("rp_start"))
            rp_l = self._safe_int(r.get("rp_len"))
            if None in (g_start, g_len, fp_s, fp_l, rp_s, rp_l):
                return 0
            overlap_fp = self._calc_overlap_interval(g_start, g_len, fp_s, fp_l)
            overlap_rp = self._calc_overlap_interval(g_start, g_len, rp_s, rp_l)
            return int(not (overlap_fp or overlap_rp))

        def _overlap_pam_row(r):
            pam_s = self._safe_int(r.get("pam_start"))
            pam_l = self._safe_int(r.get("pam_len"))
            fp_s = self._safe_int(r.get("fp_start"))
            fp_l = self._safe_int(r.get("fp_len"))
            rp_s = self._safe_int(r.get("rp_start"))
            rp_l = self._safe_int(r.get("rp_len"))
            if None in (pam_s, pam_l, fp_s, fp_l, rp_s, rp_l):
                return 0
            overlap_fp = self._calc_overlap_interval(pam_s, pam_l, fp_s, fp_l)
            overlap_rp = self._calc_overlap_interval(pam_s, pam_l, rp_s, rp_l)
            return int(not (overlap_fp or overlap_rp))

        df["overlap_protospacer"] = df.apply(_overlap_protospacer_row, axis=1)
        df["overlap_pam"] = df.apply(_overlap_pam_row, axis=1)


        # Off-target mismatch proxies:
        # Convert best percent mismatch vs background into an "importance-like" signal in [0,1]
        # (1 - mm%), without ranking/scoring aggregation.
        def _best_mm_vs_bg(query: str, bgs: list[str]) -> float:
            if not query or not bgs:
                return 100.0  # treat as max mismatch if no backgrounds given
            vals = [self._percent_mismatch(query, b[: len(query)]) for b in bgs if b]
            return min(vals) if vals else 100.0

        df["crrna_offtarget_mm_imp"] = [
            max(0.0, min(1.0, 1.0 - _best_mm_vs_bg(q, b) / 100.0))
            for q, b in zip(df["crrna"], bg_lists)
        ]
        df["fp_offtarget_mm_imp"] = [
            max(0.0, min(1.0, 1.0 - _best_mm_vs_bg(q, b) / 100.0))
            for q, b in zip(df["forward_primer"], bg_lists)
        ]
        df["rp_offtarget_mm_imp"] = [
            max(0.0, min(1.0, 1.0 - _best_mm_vs_bg(q, b) / 100.0))
            for q, b in zip(df["backward_primer"], bg_lists)
        ]

        self.df = df
        return self

    def to_csv(self, output_csv: str):
        self.df.to_csv(output_csv, index=False)
        return output_csv




