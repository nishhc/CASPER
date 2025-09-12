import numpy as np
import pandas as pd

class RankerVibe:
    FEATURE_WEIGHTS = {
        "crrna_offtarget_mm_imp": 0.08,
        "fp_offtarget_mm_imp":    0.07,
        "rp_offtarget_mm_imp":    0.07,

        "fp_3p_self_run":         0.045,
        "rp_3p_self_run":         0.045,
        "fp_rp_3p_cross_run":     0.050,

        "fp_self_dimer_pct":      0.1,
        "rp_self_dimer_pct":      0.1,
        "fp_rp_cross_dimer_pct":  0.1,
        "fp_gc_pct":              0.025,
        "rp_gc_pct":              0.025,
        "amplicon_gc_pct":        0.025,

        "seed_max_run":           0.035,
        "seed_gc_pct":            0.025,
        "guide_mfe_kcal":         0.015,
        "guide_seed_unpaired_frac": 0.03,

        "delta_tm_C":             0.08,
        "fp_tm_C":                0.03,
        "rp_tm_C":                0.03,

        "overlap_protospacer":    0.02,
        "overlap_pam":            0.02,
    }

    CATEGORY_FEATURES = {
        "Off-targets": ["crrna_offtarget_mm_imp", "fp_offtarget_mm_imp", "rp_offtarget_mm_imp"],
        "3′ end stability": ["fp_3p_self_run", "rp_3p_self_run", "fp_rp_3p_cross_run"],
        "Primer dimers / secondary structure": ["fp_self_dimer_pct", "rp_self_dimer_pct", "fp_rp_cross_dimer_pct"],
        "Guide seed & structure": ["seed_len", "seed_max_run", "seed_gc_pct", "guide_mfe_kcal", "guide_seed_unpaired_frac"],
        "Tm balance": ["delta_tm_C", "fp_tm_C", "rp_tm_C"],
        "Cross-interactions": ["overlap_protospacer", "overlap_pam"]
    }

    FEATURE_SPECS = {
        "crrna_offtarget_mm_imp": {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "fp_offtarget_mm_imp":    {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "rp_offtarget_mm_imp":    {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},

        "fp_3p_self_run":         {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "rp_3p_self_run":         {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "fp_rp_3p_cross_run":     {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},

        # H raised to 100.0 and auto_cap added so they don't all map to 0
        "fp_self_dimer_pct":      {"type": "smaller_better_lin", "L": 50, "H": 100.0, "auto_cap": 0.95},
        "rp_self_dimer_pct":      {"type": "smaller_better_lin", "L": 50, "H": 100.0, "auto_cap": 0.95},
        "fp_rp_cross_dimer_pct":  {"type": "smaller_better_lin", "L": 50, "H": 100.0, "auto_cap": 0.95},

        "fp_gc_pct":              {"type": "range_plateau", "a": 40.0, "b": 60.0, "falloff": 10.0},
        "rp_gc_pct":              {"type": "range_plateau", "a": 40.0, "b": 60.0, "falloff": 10.0},
        "amplicon_gc_pct":        {"type": "range_plateau", "a": 35.0, "b": 65.0, "falloff": 10.0},

        "seed_max_run":           {"type": "smaller_better_lin", "L": 0.0, "H": 5.0},
        "seed_gc_pct":            {"type": "range_plateau", "a": 35.0, "b": 60.0, "falloff": 10.0},
        "guide_mfe_kcal":         {"type": "bigger_better_lin", "L": -60.0, "H": 0.0},
        "guide_seed_unpaired_frac":{"type": "bigger_better_lin", "L": 0.0, "H": 1.0},

        "delta_tm_C":             {"type": "target_laplace", "t": 0.0, "s": 1.5},
        "fp_tm_C":                {"type": "range_plateau", "a": 55.0, "b": 65.0, "falloff": 7.0},
        "rp_tm_C":                {"type": "range_plateau", "a": 55.0, "b": 65.0, "falloff": 7.0},

        "overlap_protospacer":    {"type": "smaller_better_lin", "L": 0.0, "H": 1.0},
        "overlap_pam":            {"type": "smaller_better_lin", "L": 0.0, "H": 1.0},
    }

    def __init__(self, input_csv):
        self.input_csv = input_csv
        self.df = pd.read_csv(input_csv)
        wsum = sum(self.FEATURE_WEIGHTS.values())
        self.norm_weights = {k: v / wsum for k, v in self.FEATURE_WEIGHTS.items()}

    @staticmethod
    def _clip01(x):
        return np.minimum(1.0, np.maximum(0.0, x))

    def _u_smaller_better_lin(self, x, L, H):
        return self._clip01(1.0 - (x - L) / max(H - L, 1e-9))

    def _u_bigger_better_lin(self, x, L, H):
        return self._clip01((x - L) / max(H - L, 1e-9))

    def _u_smaller_better_exp(self, x, best, half_life):
        tau = max(half_life, 1e-9) / np.log(2.0)
        return self._clip01(np.exp(-(x - best) / tau))

    def _u_target_laplace(self, x, t, s):
        return np.exp(-np.abs(x - t) / max(s, 1e-9))

    def _u_range_plateau(self, x, a, b, falloff):
        inside = (x >= a) & (x <= b)
        dist = np.where(x < a, a - x, np.where(x > b, x - b, 0.0))
        u = np.where(inside, 1.0, self._clip01(1.0 - dist / max(falloff, 1e-9)))
        return u

    def _effective_bounds(self, name, x, spec):
        L = spec.get("L", None)
        H = spec.get("H", None)
        auto_q = spec.get("auto_cap", None)
        if L is not None and H is not None and auto_q:
            qH = np.nanpercentile(x, 100 * float(auto_q))
            H = max(H, qH)
            if H <= L:
                H = L + 1e-6
        return L, H

    def _feature_utility(self, name, x):
        spec = self.FEATURE_SPECS.get(name)
        if spec is None:
            q5, q95 = np.nanpercentile(x, [5, 95])
            return self._clip01((x - q5) / max(q95 - q5, 1e-9))
        t = spec["type"]
        if t == "smaller_better_lin":
            L, H = self._effective_bounds(name, x, spec)
            return self._u_smaller_better_lin(x, L, H)
        if t == "bigger_better_lin":
            L = spec["L"]; H = spec["H"]
            return self._u_bigger_better_lin(x, L, H)
        if t == "smaller_better_exp":
            return self._u_smaller_better_exp(x, spec["best"], spec["half_life"])
        if t == "target_laplace":
            return self._u_target_laplace(x, spec["t"], spec["s"])
        if t == "range_plateau":
            return self._u_range_plateau(x, spec["a"], spec["b"], spec["falloff"])
        raise ValueError(f"Unknown utility type for {name}: {t}")

    def rank(self):
        df = self.df.copy()
        for f in self.FEATURE_WEIGHTS:
            if f in df.columns:
                u = self._feature_utility(f, df[f].to_numpy(dtype=float))
                df[f"{f}__u"] = u

        df["ranking_score"] = 0.0
        for f, w in self.norm_weights.items():
            col = f"{f}__u"
            if col in df.columns:
                df["ranking_score"] += df[col] * w

        for cat, feats in self.CATEGORY_FEATURES.items():
            present = [f for f in feats if f"{f}__u" in df.columns]
            if not present:
                continue
            wsum = sum(self.norm_weights[f] for f in present if f in self.norm_weights)
            if wsum <= 0:
                continue
            cat_score = 0.0
            for f in present:
                if f in self.norm_weights:
                    cat_score += df[f"{f}__u"] * (self.norm_weights[f] / wsum)
            df[f"{cat}__subscore"] = cat_score

        df = df.sort_values("ranking_score", ascending=False).reset_index(drop=True)
        self.df = df
        return df

    def to_csv(self, output_csv):
        self.df.to_csv(output_csv, index=False)
