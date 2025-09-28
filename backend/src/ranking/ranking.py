import pandas as pd

class Ranker:
    FEATURE_WEIGHTS = {
    "crrna_offtarget_mm_imp": 0.1096,
    "fp_offtarget_mm_imp":    0.1005,
    "rp_offtarget_mm_imp":    0.1005,

    "fp_3p_self_run":         0.0593,
    "rp_3p_self_run":         0.0593,
    "fp_rp_3p_cross_run":     0.0593,

    "fp_self_dimer_pct":      0.0365,
    "rp_self_dimer_pct":      0.0365,
    "fp_rp_cross_dimer_pct":  0.0365,
    "fp_gc_pct":              0.0228,
    "rp_gc_pct":              0.0228,
    "amplicon_gc_pct":        0.0228,

    "seed_max_run":           0.0319,
    "seed_gc_pct":            0.0228,
    "guide_mfe_kcal":         0.0136,
    "guide_seed_unpaired_frac": 0.0274,

    "delta_tm_C":             0.0913,
    "fp_tm_C":                0.0274,
    "rp_tm_C":                0.0274,

    "overlap_protospacer":    0.0182,
    "overlap_pam":            0.0182,
}

    FEATURE_SPECS = {
        "crrna_offtarget_mm_imp": {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "fp_offtarget_mm_imp":    {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "rp_offtarget_mm_imp":    {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},

        "fp_3p_self_run":         {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "rp_3p_self_run":         {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},
        "fp_rp_3p_cross_run":     {"type": "smaller_better_exp", "best": 0.0, "half_life": 1.0},

        "fp_self_dimer_pct":      {"type": "smaller_better_lin", "L": 0.0, "H": 50.0},
        "rp_self_dimer_pct":      {"type": "smaller_better_lin", "L": 0.0, "H": 50.0},
        "fp_rp_cross_dimer_pct":  {"type": "smaller_better_lin", "L": 0.0, "H": 50.0},

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

    CATEGORY_FEATURES = {
    "Off-targets": ["crrna_offtarget_mm_imp", "fp_offtarget_mm_imp", "rp_offtarget_mm_imp"],
    "3′ end stability": ["fp_3p_self_run", "rp_3p_self_run", "fp_rp_3p_cross_run"],
    "Primer dimers / secondary structure": ["fp_self_dimer_pct", "rp_self_dimer_pct", "fp_rp_cross_dimer_pct"],
    "Guide seed & structure": ["seed_len", "seed_max_run", "seed_gc_pct", "guide_mfe_kcal", "guide_seed_unpaired_frac"],
    "Tm balance": ["delta_tm_C"],
    "Cross-interactions": ["overlap_protospacer", "overlap_pam"]
}

    def __init__(self, input_csv):
        self.input_csv = input_csv
        self.df = pd.read_csv(input_csv)

    def rank(self):
        df = self.df
        # Filtering
        '''df = df[
            (df["overlap_pam"] == 0) &
            (df["overlap_protospacer"] == 0) &
            (df["amplicon_len"].between(100, 220)) &
            (df["fp_len"].between(28, 36)) &
            (df["rp_len"].between(28, 36)) &
            (df["fp_gc_pct"].between(35, 65)) &
            (df["rp_gc_pct"].between(35, 65)) &
            (df["amplicon_gc_pct"].between(35, 65)) &
            (df["delta_tm_C"] <= 3.0) &
            (df["fp_3p_self_run"] <= 2) &
            (df["rp_3p_self_run"] <= 2) &
            (df["fp_gquad_3p"] == 0) &
            (df["rp_gquad_3p"] == 0) &
            (df["fp_max_hpoly_run"] <= 4) &
            (df["rp_max_hpoly_run"] <= 4) &
            (df["amplicon_max_run"] <= 5) 
            #(df["fp_tail_gc"].between(1, 3)) &
            #(df["rp_tail_gc"].between(1, 3))
        ]'''
        # Normalization
        for feature in self.FEATURE_WEIGHTS:
            if feature in df.columns:
                min_val = df[feature].min()
                max_val = df[feature].max()
                if max_val > min_val:
                    df[feature + "_norm"] = (df[feature] - min_val) / (max_val - min_val)
                else:
                    df[feature + "_norm"] = 0
        # Ranking score
        df["ranking_score"] = 0
        for feature, weight in self.FEATURE_WEIGHTS.items():
            norm_feature = feature + "_norm"
            if norm_feature in df.columns:
                df["ranking_score"] += df[norm_feature] * weight
        # Normalized score
        max_score = df["ranking_score"].max()
        if max_score > 0:
            df["ranking_score_normalized"] = df["ranking_score"]
        else:
            df["ranking_score_normalized"] = df["ranking_score"]
        df = df.sort_values(by="ranking_score_normalized", ascending=False)
        # Category subscores
        for category, features in self.CATEGORY_FEATURES.items():
            total_weight = sum(self.FEATURE_WEIGHTS[f] for f in features)
            subscore = sum(df[f] * self.FEATURE_WEIGHTS[f] for f in features if f in df.columns) / total_weight
            df[category + "_subscore"] = subscore
        return df

    def to_csv(self, output_csv):
        return self.df.to_csv(output_csv, index=False)