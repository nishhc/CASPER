import argparse
import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import classification_report
from sklearn.inspection import permutation_importance


def normalized(v):
    v = np.asarray(v, dtype=float)
    s = np.abs(v).sum()
    return (np.abs(v) / s) if s > 0 else np.zeros_like(v, dtype=float)


def softmax_temperature(p, temp=1.0):
    """Temperature-smoothed softmax on a vector of nonnegative importances."""
    p = np.asarray(p, dtype=float)

    eps = 1e-12
    logp = np.log(p + eps)
    z = (logp) / max(temp, 1e-6)
    z -= z.max()
    e = np.exp(z)
    e_sum = e.sum()
    return e / e_sum if e_sum > 0 else np.full_like(p, 1.0 / len(p))


def fit_logreg_l2(X, y, use_cv=False, random_state=42):
    if use_cv:
        clf = LogisticRegressionCV(
            Cs=20,
            cv=5,
            penalty="l2",
            solver="lbfgs",
            scoring="f1",
            n_jobs=-1,
            refit=True,
            class_weight=None,
            random_state=random_state,
            max_iter=2000
        )
    else:
        clf = LogisticRegression(
            penalty="l2",
            solver="lbfgs",
            C=1.0,
            class_weight=None,
            random_state=random_state,
            max_iter=2000
        )
    clf.fit(X, y)
    return clf


def apply_pair_smoothing(weights, names, pair_specs):

    name_to_idx = {n: i for i, n in enumerate(names)}
    w = weights.copy()
    for spec in pair_specs:
        try:
            pair, strength_str = spec.split("=")
            a, b = pair.split(":")
            alpha = float(strength_str)
            if a in name_to_idx and b in name_to_idx:
                ia, ib = name_to_idx[a], name_to_idx[b]
                m = 0.5 * (w[ia] + w[ib])
                w[ia] = (1 - alpha) * w[ia] + alpha * m
                w[ib] = (1 - alpha) * w[ib] + alpha * m
        except Exception:

            pass
    return normalized(w)


def clip_floor_cap(weights, min_w=None, max_w=None):
    w = weights.copy()
    if min_w is not None:
        w = np.maximum(w, min_w)
    if max_w is not None:
        w = np.minimum(w, max_w)
    s = w.sum()
    if s <= 0:
        w = np.full_like(w, 1.0 / len(w))
    else:
        w = w / s
    return w


def main(csv_path, use_cv=False, blend_alpha=0.5, random_state=42,
         perm_repeats=30, perm_scoring="f1",
         temp=1.0, min_w=None, max_w=None, pair_smooth=None):
    df = pd.read_csv(csv_path)

    feature_cols = [
        "fp_tm_C",
        "rp_tm_C",
        "delta_tm_C",
        "fp_self_dimer_pct",
        "rp_self_dimer_pct",
        "fp_rp_cross_dimer_pct",
        "guide_mfe_kcal",
    ]
    target_col = "Experimental_Outcome"

    missing = [c for c in feature_cols + [target_col] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing expected columns in CSV: {missing}")

    y = df[target_col].map({"Yes": 1, "No": 0})
    if y.isna().any():
        bad = df.loc[y.isna(), target_col].unique()
        raise ValueError(f"Target has unexpected values: {bad}. Expected 'Yes'/'No'.")

    X = df[feature_cols].copy()

    ok = (~X.isna().any(axis=1)) & (~y.isna())
    X, y = X.loc[ok].values, y.loc[ok].values

    clf = fit_logreg_l2(X, y, use_cv=use_cv, random_state=random_state)

    y_pred = clf.predict(X)
    print("\n=== In-sample classification report (reference only) ===")
    print(classification_report(y, y_pred, target_names=["No", "Yes"]))

    coef_importance = normalized(clf.coef_.ravel())

    r = permutation_importance(
        clf, X, y, n_repeats=perm_repeats, random_state=random_state, scoring=perm_scoring
    )
    perm_importance = normalized(r.importances_mean)

    blend_alpha = float(blend_alpha)
    blended = normalized(blend_alpha * coef_importance + (1.0 - blend_alpha) * perm_importance)

    smoothed = softmax_temperature(blended, temp=max(1.0, float(temp)))

    if pair_smooth:
        pair_list = [s for s in pair_smooth.split(",") if s.strip()]
        smoothed = apply_pair_smoothing(smoothed, feature_cols, pair_list)

    min_w = None if min_w is None else float(min_w)
    max_w = None if max_w is None else float(max_w)
    balanced = clip_floor_cap(smoothed, min_w=min_w, max_w=max_w)

    out = pd.DataFrame({
        "feature": feature_cols,
        "importance_coef_L2": coef_importance,
        "importance_perm": perm_importance,
        "importance_blended": blended,
        "importance_balanced": balanced
    }).sort_values("importance_balanced", ascending=False)

    print("\n=== Normalized & Balanced Feature Importances (sum to 1) ===")
    print(out[["feature", "importance_balanced"]])
    print("\n(Sanity check) Sum:", out["importance_balanced"].sum())

    out_path = "feature_importances_balanced_constrained_no_standardize.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved importances to: {out_path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Balanced & constrained feature importances (no scaling): L2 logistic + permutation + smoothing, normalized to sum to 1.")
    p.add_argument("--csv", required=True, help="Path to CSV with features and Experimental_Outcome.")
    p.add_argument("--cv", action="store_true", help="Use LogisticRegressionCV to tune C (L2).")
    p.add_argument("--alpha", type=float, default=0.5, help="Blend weight: 1.0=L2 coef only, 0.0=permutation only. Default 0.5.")
    p.add_argument("--perm_repeats", type=int, default=30, help="Permutation repeats (default 30).")
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument("--perm_scoring", type=str, default="f1", help="Scoring for permutation importance.")
    p.add_argument("--temp", type=float, default=2.0, help="Temperature for smoothing (>1 flattens).")
    p.add_argument("--min_w", type=float, default=0.05, help="Minimum per-feature weight after balancing.")
    p.add_argument("--max_w", type=float, default=0.30, help="Maximum per-feature weight after balancing.")
    p.add_argument("--pair_smooth", type=str, default="fp_tm_C:rp_tm_C=0.6", help="Comma-separated pair smoothing specs like 'fp_tm_C:rp_tm_C=0.6,fp_self_dimer_pct:rp_self_dimer_pct=0.3'")

    args = p.parse_args()
    main(
        csv_path=args.csv,
        use_cv=args.cv,
        blend_alpha=args.alpha,
        random_state=args.seed,
        perm_repeats=args.perm_repeats,
        perm_scoring=args.perm_scoring,
        temp=args.temp,
        min_w=args.min_w,
        max_w=args.max_w,
        pair_smooth=args.pair_smooth
    )
