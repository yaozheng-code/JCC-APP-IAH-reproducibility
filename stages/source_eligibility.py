#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Restore source-presence/NA semantics of original 25262/25263/25265 SQL.

No archived historical dataset is read. Non-recorded fluid input/urine
trajectories must remain missing (NaN), not false zero measurements.
Only new analysis-memory and output copy are changed; DB stays read-only.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import psycopg

ELIGIBILITY_SQL = r"""
WITH c AS (
  SELECT d.subject_id,d.hadm_id,d.stay_id,d.icu_intime,d.icu_outtime
  FROM mimiciv_derived.icustay_detail d
  WHERE d.stay_id = ANY(%s::bigint[])
)
SELECT c.subject_id,c.hadm_id,c.stay_id,

  EXISTS (
    SELECT 1
    FROM mimiciv_icu.inputevents p
    WHERE p.hadm_id=c.hadm_id
      AND p.amount IS NOT NULL
      AND p.amountuom='mL'
      AND p.starttime<c.icu_outtime
      AND p.endtime>c.icu_intime
  ) AS has_input_source_25262,

  EXISTS (
    SELECT 1
    FROM mimiciv_icu.outputevents p
    WHERE p.hadm_id=c.hadm_id
      AND p.value IS NOT NULL
      AND p.charttime>=c.icu_intime
      AND p.charttime<c.icu_outtime
  ) AS has_output_source_25263,

  EXISTS (
    SELECT 1
    FROM mimiciv_derived.urine_output p
    WHERE p.stay_id=c.stay_id
      AND p.urineoutput IS NOT NULL
      AND p.charttime>=c.icu_intime
      AND p.charttime<c.icu_outtime
  ) AS has_urine_source_25265

FROM c
ORDER BY c.stay_id"""

SOURCE_FIELDS = {
    "has_input_source_25262": "input",
    "has_output_source_25263": "output",
    "has_urine_source_25265": "urine",
}
WINDOWS=("pre24","0_24","24_48")

def source_flags_from_postgres(native:pd.DataFrame, *, dsn:str) -> pd.DataFrame:
    cohort=native.loc[:,["subject_id","hadm_id","stay_id"]].copy()
    for key in cohort:
        cohort[key]=pd.to_numeric(cohort[key],errors="raise").astype("int64")
    if len(cohort)!=641 or cohort.stay_id.duplicated().any():
        raise RuntimeError("The study cohort must have 641 unique ICU stays.")
    ids=[int(i) for i in cohort["stay_id"]]
    with psycopg.connect(dsn) as conn:
        with conn.cursor() as cur:
            cur.execute("SET TRANSACTION READ ONLY")
            cur.execute("EXPLAIN "+ELIGIBILITY_SQL,(ids,))
            cur.execute(ELIGIBILITY_SQL,(ids,))
            flags=pd.DataFrame(cur.fetchall(),columns=[c.name for c in cur.description])
    if len(flags)!=641 or flags.stay_id.duplicated().any():
        raise RuntimeError("Source eligibility did not return one row per ICU stay.")
    flags=flags.sort_values("stay_id").reset_index(drop=True)
    cohort=cohort.sort_values("stay_id").reset_index(drop=True)
    for key in cohort:
        if not cohort[key].equals(
            pd.to_numeric(flags[key],errors="raise").astype("int64")
        ):
            raise RuntimeError(f"Eligibility ID mismatch in {key}.")
    for col in SOURCE_FIELDS:
        if flags[col].isna().any():
            raise RuntimeError(f"Unknown source status for {col}")
        flags[col]=flags[col].astype(bool)
    return flags

def restore_source_missingness(native:pd.DataFrame, flags:pd.DataFrame):
    """Recreate source inclusion from the locally queried source tables only."""
    if not native.stay_id.equals(flags.stay_id):
        raise RuntimeError("Source flags were not aligned with native stay_id.")
    fixed=native.copy()
    report={}
    for flag,prefix in SOURCE_FIELDS.items():
        absent=~flags[flag]
        details={
            "source_eligible_stays":int(flags[flag].sum()),
            "source_ineligible_stays":int(absent.sum()),
            "newly_missing":{}
        }
        for window in WINDOWS:
            col=f"{prefix}_{window}"
            if col not in fixed:
                raise RuntimeError(f"Missing native fluid trajectory field {col}")
            was= fixed[col].isna()
            fixed.loc[absent,col]=np.nan
            details["newly_missing"][col]=int((~was&fixed[col].isna()).sum())
        report[prefix]=details
    for window in WINDOWS:
        fixed[f"net_balance_{window}"]=(
            pd.to_numeric(fixed[f"input_{window}"],errors="coerce")
            -pd.to_numeric(fixed[f"output_{window}"],errors="coerce")
        )
    fixed["net_balance_pre24_L"]=fixed["net_balance_pre24"]/1000.0
    if "positive_balance_pre24" in fixed:
        fixed["positive_balance_pre24"]=fixed.net_balance_pre24.gt(0).astype(int)
    report["net_balance_pre24_L_missing"]=int(fixed.net_balance_pre24_L.isna().sum())
    report["urine_0_24_missing"]=int(fixed.urine_0_24.isna().sum())
    # Explicitly preserve study exposure, outcomes, established clinical scores.
    invariant=("app_t0","iap_t0","mbp_pre6_last","death_28d_from_t0",
               "hospital_expire_flag","weight","aki_dynamic_first",
               "sofa_dynamic_first","vis_pre6_max")
    for col in invariant:
        if not pd.to_numeric(fixed[col],errors="coerce").equals(
            pd.to_numeric(native[col],errors="coerce")
        ):
            raise RuntimeError(f"Restoration unexpectedly altered {col}")
    return fixed, report
