#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 2: Native MIMIC-IV clinical covariates, extending previously successful
native Stage 1 results WITHOUT reading prior study extract files.

This adds first in-ICU SOFA/AKI, admission-level weight, Charlson, sepsis, invasive ventilation;
and whole-stay SOFA/AKI maximum fields for supplementary descriptions.
The original VIS is not recalculated here and is deliberately left missing.
This script is NOT a complete paper-results reproducer yet.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
DEFAULT_STAGE1=(Path.home()/"mimic-project"/"results"/"jcc_app_iah_native_stage1"
                /"native_stage1_641_clinical_partial_LOCAL_ONLY.csv")
DEFAULT_OUT=Path.home()/"mimic-project"/"results"/"jcc_app_iah_native_stage2"

STUDY_FIELDS=[
    "sofa_dynamic_first","aki_dynamic_first","charlson",
    "weight","sepsis","InvasiveVent"
]
ADDITIONAL_FIELDS=[
    "sofa_dynamic_first_time","aki_dynamic_first_time",
    "weight_first_time","weight_source_stay_id",
    "sepsis3_true_diagnostic","sofa","aki_score"
]
INT_IDS=["subject_id","hadm_id","stay_id"]
EXPECTED={"sofa_dynamic_first_sum":2774, "aki_dynamic_first_sum":116,
          "charlson_sum":2929, "weight_sum":57069.2,
          "sepsis_n":571, "InvasiveVent_n":573}

def standardized_ids(data:pd.DataFrame,context:str)->pd.DataFrame:
    data=data.copy()
    missing=[col for col in INT_IDS if col not in data]
    if missing:
        raise ValueError(f"{context} missing columns: {missing}")
    for key in INT_IDS:
        num=pd.to_numeric(data[key],errors="raise")
        if num.isna().any() or (num % 1 != 0).any():
            raise ValueError(f"{context}: invalid key {key}")
        data[key]=num.astype("int64")
    if data.stay_id.duplicated().any():
        raise ValueError(f"{context}: duplicated stay_id")
    if len(data)!=641:
        raise ValueError(f"{context}: n={len(data)}; expected 641")
    return data

def inspect_schema(cursor):
    required={
        "mimiciv_derived.icustay_detail":["subject_id","hadm_id","stay_id","icu_intime","icu_outtime"],
        "mimiciv_derived.sofa":["stay_id","starttime","endtime","sofa_24hours"],
        "mimiciv_derived.kdigo_stages":["stay_id","charttime","aki_stage"],
        "mimiciv_derived.charlson":["subject_id","hadm_id","charlson_comorbidity_index"],
        "mimiciv_derived.weight_durations":["stay_id","starttime","endtime","weight"],
        "mimiciv_derived.sepsis3":["stay_id","sepsis3"],
        "mimiciv_derived.ventilation":["stay_id","ventilation_status"],
        "mimiciv_icu.icustays":["stay_id","hadm_id"]
    }
    for table,expected in required.items():
        cursor.execute("""
          SELECT a.attname
          FROM pg_attribute a
          WHERE a.attrelid=to_regclass(%s)
          AND a.attnum>0 AND NOT a.attisdropped
        """,(table,))
        present={r[0] for r in cursor.fetchall()}
        absent=sorted(set(expected)-present)
        if absent:
            raise RuntimeError(f"Table/column missing: {table}: {absent}")
    print("Stage 2: 8 database tables/columns matched the reported schema.",flush=True)

def load_covariates(args,stage1):
    import psycopg
    sql=(HERE/"native_stage2_covariates.sql").read_text(encoding="utf-8")
    stay_ids=[int(x) for x in stage1.stay_id]
    with psycopg.connect(
        args.dsn if args.dsn else f"dbname={args.database}"
    ) as conn:
        with conn.cursor() as cur:
            cur.execute("SET TRANSACTION READ ONLY")
            cur.execute("SELECT current_database()")
            database=cur.fetchone()[0]
            inspect_schema(cur)
            cur.execute("EXPLAIN "+sql,(stay_ids,))
            print("Stage 2: SQL structure EXPLAIN PASS; extracting live local clinical fields...",flush=True)
            cur.execute(sql,(stay_ids,))
            rows=cur.fetchall()
            names=[c.name for c in cur.description]
    data=standardized_ids(pd.DataFrame(rows,columns=names),
                          "Stage 2 native SQL")
    return data,database

def combine(stage1,source):
    # Preserve all Stage1 exposures and original T0 in place; do not overwrite
    # any already calculated APP/MAP/IAP, model fields, or outcomes.
    if not set(stage1.stay_id)==set(source.stay_id):
        raise ValueError("Stage1 native cohort differs from Stage2 SQL cohort")
    left=stage1.set_index("stay_id")
    right=source.set_index("stay_id")
    for k in ("subject_id","hadm_id"):
        if not left[k].equals(right[k].reindex(left.index)):
            raise ValueError(f"Stage1/Stage2 ICU/hospital key mismatch: {k}")
    replace=STUDY_FIELDS+ADDITIONAL_FIELDS
    for col in replace:
        if col in left:
            if left[col].notna().any():
                raise ValueError(f"{col} already has values in Stage1; refusing silent overwrite")
            left=left.drop(columns=[col])
    joined=left.join(right[replace],how="inner",validate="one_to_one")
    assert len(joined)==641
    for col in STUDY_FIELDS+["sofa","aki_score","sepsis3_true_diagnostic"]:
        joined[col]=pd.to_numeric(joined[col],errors="coerce")
    return joined.reset_index().sort_values("stay_id").reset_index(drop=True)

def aggregate_qa(df):
    results={"total_icu_stays":len(df)}
    # Under the original admission-level weight rule a weight may originate
    # from another ICU stay during the same hospital admission.
    source=pd.to_numeric(df.weight_source_stay_id,errors="coerce")
    idx=pd.to_numeric(df.stay_id,errors="coerce")
    results["weight_from_other_icu_stay_count"]=int(
        (source.notna() & (source!=idx)).sum()
    )
    for field in STUDY_FIELDS:
        series=pd.to_numeric(df[field],errors="coerce")
        results[f"{field}_nonnull"]=int(series.notna().sum())
        results[f"{field}_missing"]=int(series.isna().sum())
        results[f"{field}_sum"]=float(series.sum()) if series.notna().any() else None
    for name in ("sofa","aki_score"):
        series=pd.to_numeric(df[name],errors="coerce")
        results[f"{name}_nonnull"]=int(series.notna().sum())
        results[f"{name}_sum"]=float(series.sum()) if series.notna().any() else None
    results["sepsis3_true_count_diagnostic"]=int(
        df.sepsis3_true_diagnostic.sum())
    results["vis_pre6_max_nonnull"]=int(
        pd.to_numeric(df.vis_pre6_max,errors="coerce").notna().sum()
    ) if "vis_pre6_max" in df else 0
    results["reference_aggregate_values_for_comparison_only"]=EXPECTED
    results["no_reference_patient_data_read"]=True
    return results


def optional_reference_qa(native:pd.DataFrame, csv_path:Path) -> dict:
    """Read an existing reference dataset only to count discrepancies.

    No historical reference values enter the native clinical dataset.
    No stay_id or row-level measurements are emitted to console or reports.
    """
    if not csv_path.is_file():
        return {"status":"missing_optional_reference_file"}
    ref=standardized_ids(
        pd.read_csv(csv_path, usecols=["stay_id","hadm_id","subject_id",
                                       "weight","aki_dynamic_first"]),
        "Optional reference QA dataset"
    )
    l=native.set_index("stay_id").sort_index()
    r=ref.set_index("stay_id").sort_index()
    if not l.index.equals(r.index):
        return {"status":"blocked_mismatched_reference_641_stay_ids"}
    for key in ("subject_id","hadm_id"):
        if not l[key].equals(r[key]):
            return {"status":"blocked_mismatched_reference_encounter_keys"}
    out={"status":"PASS_id_link_reference_QA_only",
         "matched_stay_ids":641,
         "comparison_alters_native_dataset":False}
    for col in ("weight","aki_dynamic_first"):
        x=pd.to_numeric(l[col],errors="coerce").astype(float)
        y=pd.to_numeric(r[col],errors="coerce").astype(float)
        both=x.notna() & y.notna()
        match=(x[both]-y[both]).abs().le(0.00001)
        out[col]={
            "compared_rows":int(both.sum()),
            "native_missing":int(x.isna().sum()),
            "reference_missing":int(y.isna().sum()),
            "different_values":int((~match).sum()),
            "either_missing":int((~both).sum()),
            "native_sum":round(float(x.sum()),2),
            "reference_sum":round(float(y.sum()),2),
            "native_minus_reference_sum":round(float(x.sum()-y.sum()),2),
        }
        if col=="aki_dynamic_first":
            differences=x[both]-y[both]
            # Aggregate difference histogram; no patient identifiers.
            counts=differences[differences.abs()>0.00001].round(3).value_counts()
            out[col]["difference_histogram"]={
                str(delta):int(count) for delta,count in counts.items()
            }
    return out

def main():
    parser=argparse.ArgumentParser(description="Local MIMIC-IV Stage2 covariate extraction")
    parser.add_argument("--stage1",type=Path,default=DEFAULT_STAGE1)
    parser.add_argument("--outdir",type=Path,default=DEFAULT_OUT)
    parser.add_argument("--database",default="mimiciv")
    parser.add_argument("--dsn",default=None)
    parser.add_argument(
        "--qc-reference-csv",type=Path,default=None,
        help="OPTIONAL historical 641-row reference dataset for aggregate-only "
             "validation; never used to build model columns."
    )
    args=parser.parse_args()
    src=args.stage1.expanduser()
    if not src.is_file():
        raise FileNotFoundError(f"Stage1 local file not found: {src}")
    print("="*65,flush=True)
    print("JCC APP/IAH | Local MIMIC-IV Stage 2 clinical covariates",flush=True)
    print("Input: your already successfully completed native Stage1 641 rows",flush=True)
    print("No historical exports. No extra cohort filtering. No database writes.",flush=True)
    print("="*65,flush=True)
    stage1=standardized_ids(pd.read_csv(src,low_memory=False),"Stage1 cohort")
    assert "t0_time" in stage1 and "app_t0" in stage1
    source,db=load_covariates(args,stage1)
    merged=combine(stage1,source)
    qc = (
        optional_reference_qa(merged,args.qc_reference_csv.expanduser())
        if args.qc_reference_csv is not None else
        {"status":"not_requested_native_only"}
    )
    report={
        "optional_reference_qa":qc,
        "weight_join":"earliest weight within same hospital admission (hadm_id)",
        "execution_utc":datetime.now(timezone.utc).isoformat(),
        "database":db,
        "source":"MIMIC-IV v3.1 local PostgreSQL (native Stage1 + Stage2)",
        "cohort":len(merged),
        "native_covar_qa":aggregate_qa(merged),
        "data_source_exports_used":False,
        "VIS_original_model_field_verified":False,
        "statistical_models_rerun":False,
        "definitions":"First SOFA/AKI values are the earliest non-null values during ICU stay; not at IAH T0. Sepsis flag follows the original study SQL row-existence rule.",
    }
    out=args.outdir.expanduser()
    out.mkdir(parents=True,exist_ok=True)
    path=out/"native_stage2_641_clinical_partial_LOCAL_ONLY.csv"
    merged.to_csv(path,index=False)
    (out/"native_stage2_summary.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",
        encoding="utf-8"
    )
    print("\n"+"="*65,flush=True)
    print("Stage 2 (local MIMIC-IV) completed",flush=True)
    print("Matched ICU stays:",len(merged),flush=True)
    for key in STUDY_FIELDS:
        vals=pd.to_numeric(merged[key],errors="coerce")
        print(f"{key}: nonnull={vals.count()}; missing={vals.isna().sum()}; "
              f"sum={vals.sum():.2f}",flush=True)
    print("Diagnostic sepsis3 TRUE:",
          report["native_covar_qa"]["sepsis3_true_count_diagnostic"],flush=True)
    print("Other-ICU weights recovered within same admission:",
          report["native_covar_qa"]["weight_from_other_icu_stay_count"],flush=True)
    print("Optional reference check (NOT used in data construction):",
          qc.get("status"),flush=True)
    if qc.get("status") == "PASS_id_link_reference_QA_only":
        print("Reference QA weight:",qc["weight"],flush=True)
        print("Reference QA first AKI:",qc["aki_dynamic_first"],flush=True)
    print("Native VIS (model score) remains unverified and missing.",flush=True)
    print("Output:",out,flush=True)
    print("Do not share any patient-level CSV.",flush=True)
    print("="*65,flush=True)

if __name__=="__main__":
    main()
