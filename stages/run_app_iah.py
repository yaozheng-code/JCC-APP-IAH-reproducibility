#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JCC APP/IAH study — analysis-level reproduction from the 641-row derived dataset.

Scope: reproduces study analyses from an existing T0-anchored patient-level
analysis dataset. It does NOT recreate that dataset directly from native
MIMIC-IV tables. It does not use albumin/furosemide treatment groups.

Requirements: see requirements.txt
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import warnings
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, roc_curve

BASE_COVARIATES = [
    "admission_age", "female", "net_balance_pre24_L",
    "albumin_pre24_last", "vis_pre6_max", "lactate_pre6_last",
    "creat_pre24_last", "k_pre24_last", "sofa_dynamic_first",
    "aki_dynamic_first", "charlson", "weight"
]
OUTCOMES = {
    "28-day mortality": "death_28d_from_t0",
    "Hospital mortality": "hospital_expire_flag"
}
EXPOSURES = {
    "Per 10 mmHg lower APP": "lower_app_per10",
    "APP <60 mmHg": "app_low60"
}
SUPPLEMENT_EXPOSURES = {
    "Per 10 mmHg lower minimum APP": "lower_app_min_per10",
    "Minimum APP <60 mmHg": "app_min_low60"
}
QUARTILES = ["Q1 lowest", "Q2", "Q3", "Q4 highest"]


def load_input(csv_path: str | None, zip_path: str | None) -> tuple[pd.DataFrame, str, str]:
    if bool(csv_path) == bool(zip_path):
        raise ValueError("Specify exactly one of --data-csv and --data-zip.")
    if csv_path:
        data = Path(csv_path).read_bytes()
        origin = str(Path(csv_path).resolve())
    else:
        with zipfile.ZipFile(zip_path) as archive:
            matches = [name for name in archive.namelist()
                       if name.endswith("data/derived/APP_IAH_analysis_dataset.csv")]
            if len(matches) != 1:
                raise ValueError("Expected exactly one data/derived/APP_IAH_analysis_dataset.csv "
                                 f"inside ZIP; found {len(matches)}.")
            data = archive.read(matches[0])
            origin = f"{Path(zip_path).resolve()}::{matches[0]}"
    frame = pd.read_csv(io.BytesIO(data))
    return frame, hashlib.sha256(data).hexdigest(), origin


def prepare(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "female" not in df and "gender" in df:
        df["female"] = df["gender"].astype(str).str.upper().str.startswith("F").astype(int)
    if "net_balance_pre24_L" not in df and "net_balance_pre24" in df:
        df["net_balance_pre24_L"] = pd.to_numeric(
            df["net_balance_pre24"], errors="coerce") / 1000
    if "app_t0" not in df and {"mbp_pre6_last", "iap_t0"}.issubset(df):
        df["app_t0"] = df["mbp_pre6_last"] - df["iap_t0"]
    if "lower_app_per10" not in df:
        df["lower_app_per10"] = -df["app_t0"] / 10
    if "app_low60" not in df:
        df["app_low60"] = (df["app_t0"] < 60).astype(int)
    if "app_low50" not in df:
        df["app_low50"] = (df["app_t0"] < 50).astype(int)
    if "app_min_pre6" not in df and {"mbp_pre6_min", "iap_t0"}.issubset(df):
        df["app_min_pre6"] = df["mbp_pre6_min"] - df["iap_t0"]
    if "lower_app_min_per10" not in df:
        df["lower_app_min_per10"] = -df["app_min_pre6"] / 10
    if "app_min_low60" not in df:
        df["app_min_low60"] = (df["app_min_pre6"] < 60).astype(int)

    required = (
        list(BASE_COVARIATES) +
        ["app_t0", "iap_t0", "mbp_pre6_last", "app_min_pre6",
         "app_low60", "app_low50", "lower_app_per10", "app_min_low60",
         "lower_app_min_per10", "death_28d_from_t0", "hospital_expire_flag",
         "app_q", "net_balance_pre24", "iap_gap_hr", "sepsis", "InvasiveVent",
         "urine_0_24", "creat_0_24_max"]
    )
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    for c in required:
        if c != "app_q":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if any(df[y].isna().any() for y in OUTCOMES.values()):
        raise ValueError("Both mortality outcome columns must be fully observed.")
    if df["app_t0"].isna().any():
        raise ValueError("APP at T0 must be available for all analysis rows.")
    if df["mbp_pre6_last"].isna().any() or df["iap_t0"].isna().any():
        raise ValueError("MBP / IAP at T0 required for AUC and APP checks.")
    for c in ["app_low60", "app_low50", "app_min_low60"]:
        if df[c].isna().any():
            raise ValueError(f"Unexpected missing exposure: {c}")
    if len(df) != 641:
        warnings.warn(f"Observed n={len(df)}; published analytic cohort has 641 rows.")
    if df["app_q"].isna().any():
        raise ValueError("app_q must be supplied with the original quartile assignment.")
    return df


def make_design(data: pd.DataFrame, exposure: str | None,
                adjusted: bool = True, component: str | None = None,
                quadratic: bool = False) -> pd.DataFrame:
    """Within-sample median imputation and a missingness flag per incomplete covariate."""
    X = pd.DataFrame(index=data.index)
    if exposure:
        X[exposure] = pd.to_numeric(data[exposure], errors="coerce")
        if X[exposure].isna().any():
            raise ValueError(f"Exposure {exposure} has missing values.")
        if quadratic:
            X["app_quadratic_term"] = X[exposure] ** 2
    if adjusted:
        fields = list(BASE_COVARIATES)
        if component:
            fields.append(component)
        for field in fields:
            values = pd.to_numeric(data[field], errors="coerce")
            if values.isna().any():
                X[field + "_missing"] = values.isna().astype(int)
                med = values.median()
                values = values.fillna(0.0 if pd.isna(med) else med)
            X[field] = values
    return sm.add_constant(X.astype(float), has_constant="add")


def fit(data: pd.DataFrame, outcome: str, exposure: str,
        adjusted: bool = True, component: str | None = "iap_t0",
        quadratic: bool = False):
    X = make_design(data, exposure, adjusted=adjusted,
                    component=component, quadratic=quadratic)
    y = data[outcome].astype(float)
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter("always")
        result = sm.Logit(y, X).fit(disp=False, maxiter=500)
    messages = list(dict.fromkeys(str(x.message) for x in observed))
    result.convergence_messages = messages
    return result, X


def effect(result, exposure: str) -> dict:
    conf = result.conf_int().loc[exposure]
    converged = bool(result.mle_retvals.get("converged", False))
    return {
        "OR": float(np.exp(result.params[exposure])),
        "CI_lower": float(np.exp(conf.iloc[0])),
        "CI_upper": float(np.exp(conf.iloc[1])),
        "p_value": float(result.pvalues[exposure]),
        "converged": converged,
        "warnings": " | ".join(result.convergence_messages)
    }


def output_effect(outcome_name, exposure_name, model_name, data,
                  outcome, exposure, adjusted, component):
    result, _ = fit(data, outcome, exposure,
                    adjusted=adjusted, component=component)
    return dict(
        outcome=outcome_name, exposure=exposure_name, model=model_name,
        n=len(data), events=int(data[outcome].sum()),
        **effect(result, exposure)
    )


def table1_quartiles(df: pd.DataFrame) -> pd.DataFrame:
    fields = [
        ("APP", "app_t0", 1), ("IAP at T0", "iap_t0", 1),
        ("MBP before T0", "mbp_pre6_last", 1),
        ("Urine 0–24 h (mL)", "urine_0_24", 0),
        ("Max creatinine 0–24 h (mg/dL)", "creat_0_24_max", 1),
        ("Pre-T0 lactate", "lactate_pre6_last", 1),
        ("First SOFA", "sofa_dynamic_first", 1),
        ("Max VIS pre6h", "vis_pre6_max", 1)
    ]
    rows = []
    for label in QUARTILES:
        d = df.loc[df["app_q"] == label]
        if len(d) == 0:
            raise ValueError(f"Expected original APP quartile {label} absent.")
        rec = {"APP quartile": label, "n": len(d)}
        for display, col, dec in fields:
            v = d[col].dropna()
            if v.empty:
                rec[display + ", median [IQR]"] = "NA"
            else:
                fmt = f"{{:.{dec}f}}"
                rec[display + ", median [IQR]"] = (
                    f"{fmt.format(v.median())} "
                    f"[{fmt.format(v.quantile(.25))}, {fmt.format(v.quantile(.75))}]"
                )
        for display, y in OUTCOMES.items():
            events = int(d[y].sum())
            rec[display] = f"{events}/{len(d)} ({100 * events/len(d):.1f}%)"
        rows.append(rec)
    return pd.DataFrame(rows)


def main_models(df):
    rows = []
    specs = [
        ("Unadjusted", False, None),
        ("Adjusted, no IAP/MBP components", True, None),
        ("Adjusted + IAP", True, "iap_t0"),
        ("Adjusted + MBP", True, "mbp_pre6_last")
    ]
    for outcome_name, y in OUTCOMES.items():
        for label, variable in EXPOSURES.items():
            for model_name, adjusted, comp in specs:
                rows.append(output_effect(
                    outcome_name, label, model_name, df, y, variable, adjusted, comp
                ))
    return pd.DataFrame(rows)


def sensitivity(df):
    subsets = [
        ("Full analytic cohort", pd.Series(True, index=df.index)),
        ("Positive pre-T0 fluid balance", df["net_balance_pre24"] > 0),
        ("IAP confirmation interval <=6 h", df["iap_gap_hr"] <= 6),
        ("Sepsis", df["sepsis"] == 1),
        ("Invasive ventilation", df["InvasiveVent"] == 1)
    ]
    rows = []
    for name, mask in subsets:
        d = df.loc[mask].copy()
        for label, ex in EXPOSURES.items():
            rows.append(output_effect(
                "28-day mortality", label, name, d,
                "death_28d_from_t0", ex, True, "iap_t0"
            ))
    return pd.DataFrame(rows)


def alternative_app(df):
    rows = []
    for outcome_name, y in OUTCOMES.items():
        for label, ex in SUPPLEMENT_EXPOSURES.items():
            rows.append(output_effect(
                outcome_name, label, "Adjusted + IAP", df, y, ex, True, "iap_t0"
            ))
    return pd.DataFrame(rows)


def threshold_50(df):
    rows = []
    for outcome_name, y in OUTCOMES.items():
        rows.append(output_effect(
            outcome_name, "APP <50 mmHg (exploratory)",
            "Adjusted + IAP", df, y, "app_low50", True, "iap_t0"
        ))
    return pd.DataFrame(rows)


def missing_data(df):
    labels = [
        "Age", "Sex", "Pre-T0 net fluid balance", "Serum albumin",
        "Maximum VIS within 6 h before T0", "Lactate", "Serum creatinine",
        "Potassium", "First dynamic SOFA score", "First dynamic AKI score",
        "Charlson comorbidity index", "Body weight", "IAP at T0"
    ]
    cols = [
        "admission_age", "female", "net_balance_pre24_L", "albumin_pre24_last",
        "vis_pre6_max", "lactate_pre6_last", "creat_pre24_last",
        "k_pre24_last", "sofa_dynamic_first", "aki_dynamic_first",
        "charlson", "weight", "iap_t0"
    ]
    records = [
        {"covariate": label, "column": col,
         "missing_n": int(df[col].isna().sum()),
         "missing_percent": 100 * df[col].isna().sum() / len(df)}
        for label, col in zip(labels, cols)
    ]
    return pd.DataFrame(records)


def quadratic_curve(df):
    result, X = fit(
        df, "death_28d_from_t0", "lower_app_per10",
        adjusted=True, component="iap_t0", quadratic=True
    )
    grid = np.linspace(
        float(df["app_t0"].quantile(.01)),
        float(df["app_t0"].quantile(.99)),
        100
    )
    risks = []
    for v in grid:
        pred_X = X.copy()
        pred_X["lower_app_per10"] = -v / 10
        pred_X["app_quadratic_term"] = (-v / 10) ** 2
        risks.append(float(result.predict(pred_X).mean()))
    curve = pd.DataFrame({
        "APP (mmHg)": grid,
        "Adjusted marginal 28-day mortality risk": risks
    })
    meta = {
        "quadratic_coefficient": float(result.params["app_quadratic_term"]),
        "quadratic_p_value_wald": float(result.pvalues["app_quadratic_term"]),
        "converged": bool(result.mle_retvals.get("converged", False)),
        "warnings": result.convergence_messages,
        "curve_definition": "Marginal standardization: mean prediction over "
                            "observed covariates, APP set to each grid value."
    }
    # The manuscript uses an APP grid spanning its 1st-99th percentiles.
    return curve, meta


def score_predictors(df):
    y = df["death_28d_from_t0"].astype(int).to_numpy()
    clinical_X = make_design(df, exposure=None, adjusted=True, component=None)
    clinical_fit = sm.Logit(y, clinical_X).fit(disp=False, maxiter=500)
    plus_X = clinical_X.copy()
    plus_X["app_t0"] = df["app_t0"].astype(float)
    plus_fit = sm.Logit(y, plus_X).fit(disp=False, maxiter=500)
    scores = {
        "IAP alone": df["iap_t0"].to_numpy(dtype=float),
        "MBP alone": -df["mbp_pre6_last"].to_numpy(dtype=float),
        "APP alone": -df["app_t0"].to_numpy(dtype=float),
        "Clinical covariates": np.asarray(clinical_fit.predict(clinical_X)),
        "Clinical covariates + APP": np.asarray(plus_fit.predict(plus_X)),
    }
    return y, scores


def auc_analysis(df, bootstrap, seed):
    y, scores = score_predictors(df)
    names = list(scores)
    estimates = {name: roc_auc_score(y, scores[name]) for name in names}
    rng = np.random.default_rng(seed)
    draws = {name: [] for name in names}
    diffs = {
        "APP minus IAP": ("APP alone", "IAP alone"),
        "APP minus MBP": ("APP alone", "MBP alone"),
        "Clinical+APP minus Clinical": (
            "Clinical covariates + APP", "Clinical covariates"
        )
    }
    drawdiffs = {name: [] for name in diffs}
    n = len(y)
    for _ in range(bootstrap):
        sample = rng.integers(0, n, size=n)
        by = y[sample]
        if by.min() == by.max():
            continue
        aucs = {name: roc_auc_score(by, scores[name][sample]) for name in names}
        for name, val in aucs.items():
            draws[name].append(val)
        for name, (left, right) in diffs.items():
            drawdiffs[name].append(aucs[left] - aucs[right])
    rows = []
    for name in names:
        lo, hi = np.percentile(draws[name], [2.5, 97.5])
        rows.append({"metric": name, "estimate": estimates[name],
                     "bootstrap_CI_low": lo, "bootstrap_CI_high": hi,
                     "number_bootstraps": len(draws[name])})
    for name, (left, right) in diffs.items():
        lo, hi = np.percentile(drawdiffs[name], [2.5, 97.5])
        rows.append({
            "metric": name, "estimate": estimates[left] - estimates[right],
            "bootstrap_CI_low": lo, "bootstrap_CI_high": hi,
            "number_bootstraps": len(drawdiffs[name])
        })
    # Bootstrap here resamples patient-level fixed fitted scores. It does not
    # refit the clinical regression at each replicate and is not external validation.
    return pd.DataFrame(rows), (y, scores)


def make_figures(df, quartiles, curve, auc_scores, path):
    path.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7.0, 4.3))
    plt.hist(df["app_t0"].dropna(), bins=28, edgecolor="white")
    plt.axvline(60, linestyle="--", label="APP = 60 mmHg")
    plt.xlabel("APP at T0 (mmHg)")
    plt.ylabel("Patients (n)")
    plt.title("Figure 1. Distribution of APP")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path/"figure1_APP_distribution.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7.0, 4.3))
    x = np.arange(len(QUARTILES))
    q_d = [df.loc[df.app_q.eq(label)] for label in QUARTILES]
    vals = [100*d.death_28d_from_t0.mean() for d in q_d]
    hos = [100*d.hospital_expire_flag.mean() for d in q_d]
    width = .37
    plt.bar(x-width/2, vals, width=width, label="28-day mortality")
    plt.bar(x+width/2, hos, width=width, label="Hospital mortality")
    plt.xticks(x, QUARTILES)
    plt.ylabel("Mortality (%)")
    plt.title("Figure 2. Mortality by APP quartile")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path/"figure2_mortality_by_quartile.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7.0, 4.3))
    plt.plot(curve["APP (mmHg)"], curve.iloc[:, 1])
    plt.axvline(60, linestyle="--", label="APP = 60 mmHg")
    plt.xlabel("APP at T0 (mmHg)")
    plt.ylabel("Adjusted 28-day mortality risk")
    plt.title("Figure 3. Adjusted APP–mortality association")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path/"figure3_adjusted_curve.png", dpi=180)
    plt.close()

    y, scores = auc_scores
    plt.figure(figsize=(6.6, 5.3))
    for name in ["APP alone", "IAP alone", "MBP alone"]:
        fpr, tpr, _ = roc_curve(y, scores[name])
        plt.plot(fpr, tpr, label=f"{name} (AUC {roc_auc_score(y,scores[name]):.3f})")
    plt.plot([0, 1], [0, 1], linestyle=":", label="Chance")
    plt.xlabel("1 − Specificity")
    plt.ylabel("Sensitivity")
    plt.title("Figure 4. Mortality discrimination")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(path/"figure4_ROC.png", dpi=180)
    plt.close()


def report_and_compare(df, main, alt, threshold, quadratic, auc,
                       sha256, origin, outfile):
    # Published benchmark values are aggregate-only and included in the code package.
    expected = {
        "n": 641, "death28": 225, "hospital_deaths": 225,
        "app_lt60": 367, "app_lt50": 180,
        "primary_or_28day": 2.03,
        "primary_or_hospital": 2.00,
        "continuous_or_28day": 1.22,
        "exploratory50_or_28day": 1.32,
        "quadratic_p": .774,
        "auc_app": .653,
        "auc_iap": .527,
        "auc_mbp": .649
    }
    sub = main
    a = sub.loc[(sub.outcome=="28-day mortality")&
                (sub.exposure=="APP <60 mmHg")&
                (sub.model=="Adjusted + IAP"),"OR"].iloc[0]
    b = sub.loc[(sub.outcome=="Hospital mortality")&
                (sub.exposure=="APP <60 mmHg")&
                (sub.model=="Adjusted + IAP"),"OR"].iloc[0]
    c = sub.loc[(sub.outcome=="28-day mortality")&
                (sub.exposure=="Per 10 mmHg lower APP")&
                (sub.model=="Adjusted + IAP"),"OR"].iloc[0]
    low50 = threshold.loc[
        threshold.outcome=="28-day mortality","OR"].iloc[0]
    point_auc = dict(zip(auc["metric"], auc["estimate"]))
    actual = {
        "n": int(len(df)),
        "death28": int(df.death_28d_from_t0.sum()),
        "hospital_deaths": int(df.hospital_expire_flag.sum()),
        "app_lt60": int(df.app_low60.sum()),
        "app_lt50": int(df.app_low50.sum()),
        "primary_or_28day": float(a),
        "primary_or_hospital": float(b),
        "continuous_or_28day": float(c),
        "exploratory50_or_28day": float(low50),
        "quadratic_p": float(quadratic["quadratic_p_value_wald"]),
        "auc_app": float(point_auc["APP alone"]),
        "auc_iap": float(point_auc["IAP alone"]),
        "auc_mbp": float(point_auc["MBP alone"])
    }
    rows = []
    for name, reference in expected.items():
        value = actual[name]
        tol = 0 if isinstance(reference, int) else (.003 if name.startswith("auc_") else .007)
        rows.append({"check": name, "observed": value, "published": reference,
                     "tolerance": tol, "status": "PASS" if abs(value-reference)<=tol else "REVIEW"})
    report = {
        "input_origin": origin,
        "input_csv_sha256": sha256,
        "scope": "Analysis-level reproduction from an already-derived patient dataset",
        "study": "JCC APP at confirmed IAH, MIMIC-IV v3.1",
        "n": len(df),
        "outcomes": actual,
        "quadratic_model": quadratic,
        "audit_statuses": rows,
        "important_note": (
            "The original 641-row dataset must be reconstructed independently "
            "for database-to-results reproducibility. Confidence intervals from a "
            "freshly seeded bootstrap can differ slightly from the submitted table."
        )
    }
    (outfile/"report_summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(outfile/"benchmark_checks.csv", index=False)
    return rows


def save_tables(table_dir, mapping):
    table_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in mapping.items():
        frame.to_csv(table_dir/name, index=False, encoding="utf-8-sig")


def main():
    parser = argparse.ArgumentParser(
        description="Reproduce JCC T0-anchored APP/IAH statistics from the derived dataset."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--data-csv", type=str, help="Path to APP_IAH_analysis_dataset.csv")
    group.add_argument("--data-zip", type=str, help="Path to existing full reproducibility ZIP")
    parser.add_argument("--outdir", type=str, default="outputs/reproduced")
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.bootstrap < 50:
        parser.error("--bootstrap must be at least 50.")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df, sha256, origin = load_input(args.data_csv, args.data_zip)
    df = prepare(df)

    # Verify the APP identity using the supplied T0 measurements.
    app_diff = (df["mbp_pre6_last"]-df["iap_t0"]-df["app_t0"]).abs()
    if app_diff.max() > 1e-6:
        raise ValueError(f"APP ≠ MBP − IAP in input; maximum difference {app_diff.max()}.")

    quartiles = table1_quartiles(df)
    models = main_models(df)
    t2 = models.loc[models["model"]=="Adjusted + IAP"].copy()
    t3 = sensitivity(df)
    s2 = alternative_app(df)
    threshold = threshold_50(df)
    missing = missing_data(df)
    curve, quad_info = quadratic_curve(df)
    auc_table, auc_scores = auc_analysis(df, args.bootstrap, args.seed)
    save_tables(out/"tables", {
        "Table1_quartiles.csv": quartiles,
        "Table2_main_adjusted.csv": t2,
        "Table3_sensitivity.csv": t3,
        "Supplementary_S1_all_models.csv": models,
        "Supplementary_S2_alternative_APP.csv": s2,
        "Supplementary_S3_AUC_bootstrap.csv": auc_table,
        "Supplementary_S4_adjusted_curve_100_points.csv": curve,
        "Supplementary_S5_missing_data.csv": missing,
        "Exploratory_APP_below50.csv": threshold
    })
    make_figures(df, quartiles, curve, auc_scores, out/"figures")
    checks = report_and_compare(
        df, models, s2, threshold, quad_info, auc_table, sha256, origin, out
    )
    counts = {st: sum(r["status"]==st for r in checks) for st in ("PASS","REVIEW")}
    print(f"Analysis completed: N={len(df)}, 28-day deaths="
          f"{int(df.death_28d_from_t0.sum())}, hospital deaths="
          f"{int(df.hospital_expire_flag.sum())}.")
    print(f"Primary 28-day adjusted OR for APP <60 mmHg: "
          f"{t2.loc[(t2.outcome=='28-day mortality') & (t2.exposure=='APP <60 mmHg'),'OR'].iloc[0]:.4f}")
    print(f"Quadratic APP term P value: {quad_info['quadratic_p_value_wald']:.4f}")
    print(f"Benchmark: {counts['PASS']} PASS; {counts['REVIEW']} REVIEW.")
    print(f"See generated outputs: {out.resolve()}")
    if counts["REVIEW"]:
        print("Please inspect benchmark_checks.csv; no unverified claim of exact replication is made.")


if __name__ == "__main__":
    main()
