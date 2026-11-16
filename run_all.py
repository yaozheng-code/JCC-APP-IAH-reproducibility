#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""JCC APP–IAH | One-command native MIMIC-IV v3.1 approximate reproduction.

1) Rebuild 819 candidate IAH ICU stays and 641 first adult ICU stays in
   read-only local PostgreSQL (exact cohort SQL from the validated notebook).
2) Extract 10 native source trajectories, reconstruct T0/APP/outcomes.
3) Extract 6 clinical covariates directly from MIMIC derived tables.
4) Recreate six-drug VIS using MIMIC native dose rates and patient weight.
5) Restore input/output/urine source-existence NA semantics from native SQL.
6) Call original UNMODIFIED manuscript analysis implementation, producing
   tables, curves, figures, and 13 AGGREGATE benchmark comparisons.

No old 14 source exports and no archived 641-row analysis dataset are needed.
NO historical data values are used to fill any new clinical field.
No database writes. All patient-level data remain on the user's local machine.

This pipeline reports aggregate model results as observed; numerical variation
between independent derived-concept installations is possible. The aggregate
checks are not a test of exact patient-level equality.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parent
STAGES=ROOT/"stages"
sys.path.insert(0,str(STAGES))
import numpy as np
import pandas as pd
import psycopg
from compute_vis import candidate_hourly_drug_scores
from source_eligibility import source_flags_from_postgres, restore_source_missingness

DEFAULT_OUTPUT=Path.home()/"mimic-project"/"results"/"jcc_app_iah_native_approx_reproduction"

def safe_ids(frame:pd.DataFrame,name:str) -> pd.DataFrame:
    frame=frame.copy()
    for key in ("subject_id","hadm_id","stay_id"):
        if key not in frame:
            raise RuntimeError(f"{name} missing expected ID column {key}")
        x=pd.to_numeric(frame[key],errors="raise")
        if x.isna().any() or (x%1!=0).any():
            raise RuntimeError(f"{name} contains invalid {key}")
        frame[key]=x.astype("int64")
    if frame.stay_id.duplicated().any():
        raise RuntimeError(f"{name} includes duplicate stay_id values")
    return frame.sort_values("stay_id").reset_index(drop=True)

def execute_python(script:Path, args:list[str]):
    cmd=[sys.executable,str(script)]+args
    print("\n"+"="*72,flush=True)
    print(f"RUN: {script.name}",flush=True)
    print("="*72,flush=True)
    subprocess.run(cmd,check=True)

def build_cohort(*, dsn:str, output:Path):
    sql819=(STAGES/"01_cohort_819_candidates.sql").read_text(encoding="utf-8")
    sql641=(STAGES/"02_cohort_641_first_adult_stay.sql").read_text(encoding="utf-8")
    with psycopg.connect(dsn) as conn:
        with conn.cursor() as cur:
            cur.execute("SET TRANSACTION READ ONLY")
            cur.execute("SELECT current_database()")
            database=cur.fetchone()[0]
            print("Database:",database,flush=True)
            cur.execute("EXPLAIN "+sql819)
            print("Candidate cohort SQL: EXPLAIN PASS.",flush=True)
            cur.execute(sql819)
            a=pd.DataFrame(cur.fetchall(),columns=[c.name for c in cur.description])
            if not len(a):
                raise RuntimeError("No ICU stays in the IAH candidate cohort")
            a=safe_ids(a,"IAH candidates")
            candidate_ids=a.stay_id.tolist()
            cur.execute("EXPLAIN "+sql641,(candidate_ids,))
            cur.execute(sql641,(candidate_ids,))
            b=pd.DataFrame(cur.fetchall(),columns=[c.name for c in cur.description])
            b=safe_ids(b,"final eligible ICU stays")
    if not set(b.stay_id).issubset(set(a.stay_id)):
        raise RuntimeError("The final cohort is not a subset of the candidate cohort")
    print(f"Source-derived IAH candidates: {len(a)} (reference: 819)",flush=True)
    print(f"Source-derived first-adult ICU stays: {len(b)} (reference: 641)",flush=True)
    output.mkdir(parents=True,exist_ok=True)
    a.to_csv(output/"IAH_candidates_LOCAL_ONLY.csv",index=False)
    b_path=output/"Final_cohort_641_LOCAL_ONLY.csv"
    b.to_csv(b_path,index=False)
    provenance={
        "source":"MIMIC-IV v3.1 local PostgreSQL",
        "database":database,
        "candidate_n":int(len(a)),
        "final_n":int(len(b)),
        "candidate_count_matches_reference_819":bool(len(a)==819),
        "final_count_matches_reference_641":bool(len(b)==641),
        "two_step_cohort_sql_sha256":{
            "candidate":hashlib.sha256(sql819.encode()).hexdigest(),
            "first_stays":hashlib.sha256(sql641.encode()).hexdigest()
        },
        "uses_prior_cohort_patient_data":False,
        "read_only_postgresql":True,
    }
    (output/"cohort_summary.json").write_text(
        json.dumps(provenance,indent=2,ensure_ascii=False)+"\n",
        encoding="utf-8"
    )
    if len(b)!=641:
        raise RuntimeError(
            "The independent local cohort is NOT 641. Stop rather than "
            "force the published cohort size. Check MIMIC-IV version, "
            "IAP itemid 223763, derived schema, and first-encounter rules. "
            f"Actual cohort has {len(b)} ICU stays."
        )
    if len(a)!=819:
        print("REVIEW: The 819-candidate benchmark differs; inspect local setup.",
              flush=True)
    return b_path,provenance

def integrate_native(*, stage2_path:Path,drugfile:Path,dsn:str,
                     destination:Path):
    from native_stage2 import standardized_ids
    native=standardized_ids(
        pd.read_csv(stage2_path,low_memory=False),
        "Local Stage2 641 ICU stays"
    )
    if native["vis_pre6_max"].notna().any():
        raise RuntimeError(
            "Native Stage2 VIS already populated; refusing untracked overwrite."
        )
    drugs=pd.read_csv(drugfile,low_memory=False)
    candidate,vis_summary=candidate_hourly_drug_scores(native,drugs)
    candidate=standardized_ids(candidate,"Native VIS candidate source")
    native=native.sort_values("stay_id").reset_index(drop=True)
    candidate=candidate.sort_values("stay_id").reset_index(drop=True)
    for k in ("subject_id","hadm_id","stay_id"):
        if not native[k].equals(candidate[k]):
            raise RuntimeError("Native VIS source patient/admission key mismatch")
    for summary in ("max","last","mean"):
        source=f"vis_6drug_pre6_{summary}_candidate"
        dest=f"vis_pre6_{summary}"
        if dest in native and native[dest].notna().any():
            raise RuntimeError(f"Native {dest} field unexpectedly already populated")
        native[dest]=pd.to_numeric(candidate[source],errors="coerce").astype("float64")
    # Preserve standard MIMIC VIS source missingness: no eligible vasoactive
    # dose anywhere in stay -> NaN; source bins of zero can be true zero.
    print("Native six-drug VIS nonmissing:",
          int(native["vis_pre6_max"].notna().sum()),flush=True)
    print("Native six-drug VIS missing:",
          int(native["vis_pre6_max"].isna().sum()),flush=True)

    flags=source_flags_from_postgres(native,dsn=dsn)
    native,source_summary=restore_source_missingness(native,flags)
    print("Native source eligibility and missingness:",source_summary,flush=True)
    if native.app_t0.isna().any():
        raise RuntimeError("Local APP exposure missing after cohort construction")
    if int(native.app_low60.sum())!=367:
        print("REVIEW: APP<60 count differs from known original reference 367.",
              flush=True)

    destination.parent.mkdir(parents=True,exist_ok=True)
    native.to_csv(destination,index=False)
    metadata={
        "native_source_only":True,
        "patient_level_model_data_path":str(destination),
        "VIS_definition":"6-drug, includes vasopressin U/h/(60*study_weight_kg)",
        "VIS_original_source_values_recovered_exactly":False,
        "VIS_dose_source_stats":vis_summary,
        "VIS_nonmissing":int(native["vis_pre6_max"].notna().sum()),
        "VIS_missing":int(native["vis_pre6_max"].isna().sum()),
        "fluid_urine_presence_checks":source_summary,
        "n":int(len(native)),
        "APP_lt60":int(native.app_low60.sum()),
        "APP_lt50":int(native.app_low50.sum()),
        "death28":int(native.death_28d_from_t0.sum()),
        "hospital_deaths":int(native.hospital_expire_flag.sum()),
    }
    (destination.parent/"native_clinical_source_summary.json").write_text(
        json.dumps(metadata,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"
    )
    return destination,metadata

def main():
    parser=argparse.ArgumentParser(
        description="Approximate end-to-end MIMIC-IV v3.1 local reproducibility, code+aggregate benchmarks only."
    )
    parser.add_argument("--database",default="mimiciv",
                        help="Local PostgreSQL database name; default mimiciv")
    parser.add_argument("--dsn",default=None,
                        help="Optional PostgreSQL DSN; recommend libpq PG* env or .pgpass, "
                             "not a password on command line")
    parser.add_argument("--outdir",type=Path,default=DEFAULT_OUTPUT,
                        help="Local path for all PRIVATE patient-level results")
    parser.add_argument("--bootstrap",type=int,default=200,
                        help="Bootstrap replicates; >=50; use 500 for original analysis default")
    parser.add_argument("--seed",type=int,default=42)
    args=parser.parse_args()
    if args.bootstrap<50:
        parser.error("--bootstrap must be at least 50")
    # Local files should not inadvertently become world-readable on Unix.
    os.umask(0o077)
    output=args.outdir.expanduser().resolve()
    output.mkdir(parents=True,exist_ok=True,mode=0o700)
    dsn=args.dsn or f"dbname={args.database}"

    print("="*76,flush=True)
    print("JCC APP–IAH | FULL LOCAL MIMIC-IV v3.1 REPRODUCIBILITY PIPELINE",flush=True)
    print("Code-only package; no historic 14 exports or prior patient CSV required.",flush=True)
    print("All PostgreSQL operations READ ONLY. No PHI is shipped.",flush=True)
    print("Status: approximate reproducibility; not exact row-level identity.",flush=True)
    print("="*76,flush=True)

    cohort_path,cohort_info=build_cohort(dsn=dsn,output=output/"01_cohort")
    stage1_dir=output/"02_native_source_trajectories"
    stage2_dir=output/"03_native_covariates"
    analysis_input_dir=output/"04_local_clinical_dataset"
    analysis_dir=output/"05_statistical_results"

    execute_python(
        STAGES/"native_stage1.py",
        ["--database",args.database,"--cohort",str(cohort_path),
         "--outdir",str(stage1_dir)]+(["--dsn",args.dsn] if args.dsn else [])
    )
    execute_python(
        STAGES/"native_stage2.py",
        ["--database",args.database,
         "--stage1",str(stage1_dir/"native_stage1_641_clinical_partial_LOCAL_ONLY.csv"),
         "--outdir",str(stage2_dir)]+(["--dsn",args.dsn] if args.dsn else [])
    )

    data_path,meta=integrate_native(
        stage2_path=stage2_dir/"native_stage2_641_clinical_partial_LOCAL_ONLY.csv",
        drugfile=stage1_dir/"native_source_trajectories"/"native_25267_norepinephrine_mcgkgmin.csv.gz",
        dsn=dsn,
        destination=analysis_input_dir/"APP_IAH_native_analysis_641_LOCAL_ONLY.csv"
    )

    execute_python(
        STAGES/"run_app_iah.py",
        ["--data-csv",str(data_path),"--outdir",str(analysis_dir),
         "--bootstrap",str(args.bootstrap),"--seed",str(args.seed)]
    )

    benchmark_path=analysis_dir/"benchmark_checks.csv"
    checks=pd.read_csv(benchmark_path)
    counts={
        "PASS":int(checks.status.eq("PASS").sum()),
        "REVIEW":int(checks.status.eq("REVIEW").sum()),
    }
    summary=json.loads((analysis_dir/"report_summary.json").read_text(encoding="utf8"))
    final={
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "study":"APP at confirmed IAH and mortality in critically ill patients",
        "source":"Local authorized MIMIC-IV v3.1 PostgreSQL, not an external data platform",
        "stage_count":5,
        "cohort":cohort_info,
        "native_clinical_summary":meta,
        "analysis_13_aggregate_benchmarks":counts,
        "analysis_benchmark_rows":checks.to_dict(orient="records"),
        "model_aggregate_estimates":summary.get("outcomes"),
        "bootstrap_repetitions":args.bootstrap,
        "bootstrap_seed":args.seed,
        "source_patient_data_in_distributed_ZIP":False,
        "historical_14_source_exports_used":False,
        "historical_641_analysis_file_used":False,
        "clinical_comparison_scope":"Aggregate numerical reproducibility; individual measurements may vary between implementations.",
        "appropriate_claim":"Approximate independent reproducibility, not identical patient-level output",
        "warning":"PASS/REVIEW is an aggregate numerical comparison; do not infer exact row-level equality.",
    }
    (output/"REPRODUCTION_SUMMARY_AGGREGATE_ONLY.json").write_text(
        json.dumps(final,indent=2,ensure_ascii=False,default=str)+"\n",
        encoding="utf-8")
    print("\n"+"="*76,flush=True)
    print("FINISHED: Local MIMIC-IV database -> clinical dataset -> study statistics",flush=True)
    print("Aggregate benchmark checks:",counts,flush=True)
    print("Source-derived first-adult IAH cohort:",cohort_info["final_n"],flush=True)
    print("Native 28-day deaths:",meta["death28"],flush=True)
    print("Local output directory:",output,flush=True)
    print("SHARE ONLY aggregate tables/figures and source code, NOT patient CSVs.",flush=True)
    print("This is not a certification of 100% patient-level equivalence.",flush=True)
    print("="*76,flush=True)

if __name__=="__main__":
    main()
