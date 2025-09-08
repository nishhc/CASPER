import pandas as pd

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
dftwo = pd.read_csv("pairs_features.csv")
df = dftwo[
    (dftwo["overlap_pam"] == 0) &
    (dftwo["overlap_protospacer"] == 0) &
    (dftwo["amplicon_len"].between(100, 220)) &
    (dftwo["fp_len"].between(28, 36)) &
    (dftwo["rp_len"].between(28, 36)) &
    (dftwo["fp_gc_pct"].between(35, 65)) &
    (dftwo["rp_gc_pct"].between(35, 65)) &
    (dftwo["amplicon_gc_pct"].between(35, 65)) &
    (dftwo["delta_tm_C"] <= 3.0) &
    (dftwo["fp_3p_self_run"] <= 2) &
    (dftwo["rp_3p_self_run"] <= 2) &
    (dftwo["fp_gquad_3p"] == 0) &
    (dftwo["rp_gquad_3p"] == 0) &
    (dftwo["fp_max_hpoly_run"] <= 4) &
    (dftwo["rp_max_hpoly_run"] <= 4) &
    (dftwo["amplicon_max_run"] <= 5) &
    (dftwo["fp_tail_gc"].between(1, 3)) &
    (dftwo["rp_tail_gc"].between(1, 3))
]
for feature in FEATURE_WEIGHTS:
    if feature in df.columns:
        min_val = df[feature].min()
        max_val = df[feature].max()
        if max_val > min_val:
            df[feature + "_norm"] = (df[feature] - min_val) / (max_val - min_val)
        else:
            df[feature + "_norm"] = 0

df["ranking_score"] = 0
for feature, weight in FEATURE_WEIGHTS.items():
    norm_feature = feature + "_norm"
    if norm_feature in df.columns:
        df["ranking_score"] += df[norm_feature] * weight


max_score = df["ranking_score"].max()
if max_score > 0:
    df["ranking_score_normalized"] = df["ranking_score"]
else:
    df["ranking_score_normalized"] = df["ranking_score"]
df = df.sort_values(by="ranking_score_normalized", ascending=False)
df.to_csv("pairs_ranked.csv", index=False)


CATEGORY_FEATURES = {
    "Off-targets": ["crrna_offtarget_mm_imp", "fp_offtarget_mm_imp", "rp_offtarget_mm_imp"],
    "3′ end stability": ["fp_3p_self_run", "rp_3p_self_run", "fp_rp_3p_cross_run"],
    "Primer dimers / secondary structure": ["fp_self_dimer_pct", "rp_self_dimer_pct", "fp_rp_cross_dimer_pct"],
    "Guide seed & structure": ["seed_len", "seed_max_run", "seed_gc_pct", "guide_mfe_kcal", "guide_seed_unpaired_frac"],
    "Tm balance": ["delta_tm_C"],
    "Cross-interactions": ["overlap_protospacer", "overlap_pam"],
    "Conservation": ["guide_conservation", "fp_conservation", "rp_conservation"],
}


for category, features in CATEGORY_FEATURES.items():
    total_weight = sum(FEATURE_WEIGHTS[f] for f in features)
    subscore = sum(df[f] * FEATURE_WEIGHTS[f] for f in features if f in df.columns) / total_weight
    df[category + "_subscore"] = subscore


print(df[[cat + "_subscore" for cat in CATEGORY_FEATURES.keys()]].head())