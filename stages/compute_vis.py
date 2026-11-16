#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MIMIC-IV six-drug vasoactive-inotropic score (VIS) source reconstruction.

Six-drug study formula:
dopamine + dobutamine + 100*epinephrine + 100*norepinephrine
+ 10*milrinone + 10000*vasopressin (U/kg/min).

MIMIC standard 'vasoactive_agent' vasopressin is units/hour; divide by
60 minutes/hour and study bodyweight (kg). The seven-drug exploratory
variant including 10*phenylephrine is computed for QC but NEVER used
in the manuscript model. Different implementations of vasoactive dose
normalization may introduce small numerical variation in independent runs.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

KEYS=["subject_id","hadm_id","stay_id"]
DRUGS=[
    "dopamine_mcgkgmin",
    "dobutamine_mcgkgmin",
    "epinephrine_mcgkgmin",
    "norepinephrine_mcgkgmin",
    "milrinone_mcgkgmin",
    "phenylephrine_mcgkgmin",
    "vasopressin_units_hour",
]

def candidate_hourly_drug_scores(stage2, raw):
    """Compute 6-drug/7-drug VIS candidates; DO NOT choose an original score.

    MIMIC-IV vasoactive_agent.sql: vasopressin=units/hour, other drug
    rates=mcg/kg/min. Published VIS score definitions expect
    vasopressin in units/kg/min, so divide by (60*weight_kg).

    These two variants differ only in the addition of 10*phenylephrine.
    The actual study's custom VIS formula has NOT yet been proven.
    """
    required=KEYS+["t_hours_from_icu_intime"]+DRUGS
    missing=[col for col in required if col not in raw]
    if missing:
        raise RuntimeError(f"Native drug-dose extract missing columns: {missing}")
    raw=raw.loc[:,required].copy()
    if raw.empty:
        raise RuntimeError("Native drug-dose input is empty")
    for k in KEYS:
        raw[k]=pd.to_numeric(raw[k],errors="raise").astype("int64")
    if not set(raw.stay_id).issubset(set(stage2.stay_id)):
        raise RuntimeError("Drug-dose source contains unknown stays")
    if raw.duplicated(["stay_id","t_hours_from_icu_intime"]).any():
        raise RuntimeError("Duplicate native drug dose bin per stay/hour")
    raw["t_hours_from_icu_intime"]=pd.to_numeric(
        raw["t_hours_from_icu_intime"],errors="raise").astype("float64")
    for col in DRUGS:
        raw[col]=pd.to_numeric(raw[col],errors="coerce").astype("float64")
    if raw[DRUGS].lt(-1e-7).any().any():
        raise RuntimeError("Negative vasoactive-agent dose; review units/source")

    stage2=stage2.copy()
    stage2["weight"]=pd.to_numeric(stage2["weight"],errors="coerce").astype("float64")
    stage2["t0_hr"]=pd.to_numeric(stage2["t0_hr"],errors="coerce").astype("float64")
    if not (stage2["weight"].notna().all() and stage2["weight"].gt(0).all()):
        raise RuntimeError("Invalid/missing weight; stop rather than invent a vasopressin conversion")
    if not stage2["t0_hr"].notna().all():
        raise RuntimeError("Missing IAH T0 hour")
    mapping=stage2.set_index("stay_id")[["weight","t0_hr","hadm_id","subject_id"]]
    raw=raw.join(mapping[["weight","t0_hr"]],on="stay_id")
    if raw["weight"].isna().any() or raw["t0_hr"].isna().any():
        raise RuntimeError("Native drug source cannot link to Stage2")
    # Check medication source keys for every stay, not just matching stay_id.
    keys=raw[["subject_id","hadm_id","stay_id"]].drop_duplicates().set_index("stay_id")
    test=keys.join(mapping[["subject_id","hadm_id"]],rsuffix="_stage2",validate="one_to_one")
    if not (test["subject_id"]==test["subject_id_stage2"]).all():
        raise RuntimeError("Drug-dose source subject_id mismatch")
    if not (test["hadm_id"]==test["hadm_id_stage2"]).all():
        raise RuntimeError("Drug-dose source hadm_id mismatch")

    doses=raw[DRUGS].fillna(0.0)
    # Rate values are rounded to six decimals by the previously executed
    # Stage1 source SQL; no re-normalization of drug-specific rates occurs.
    standard_6 = (
        doses["dopamine_mcgkgmin"]
        + doses["dobutamine_mcgkgmin"]
        + 100.0*doses["epinephrine_mcgkgmin"]
        + 100.0*doses["norepinephrine_mcgkgmin"]
        + 10.0*doses["milrinone_mcgkgmin"]
        + (10000.0/60.0)*doses["vasopressin_units_hour"]/raw["weight"]
    )
    raw["vis_6drug_candidate"]=standard_6
    raw["vis_7drug_candidate"]=standard_6+10.0*doses["phenylephrine_mcgkgmin"]
    raw["drug_rate_any_positive"]=doses.gt(1e-9).any(axis=1)
    active_stays=set(raw.loc[raw.drug_rate_any_positive,"stay_id"].astype("int64"))
    raw["rel_hours"]=raw["t_hours_from_icu_intime"]-raw["t0_hr"]
    pre6=raw["rel_hours"].between(-6.0,0.0,inclusive="both")
    # Historical 25267 source eligibility required a non-null drug score
    # somewhere in the ICU stay. This native approximation uses positive
    # drug-dose presence somewhere; the discrepancy is explicitly audited.
    relevant=raw.loc[pre6 & raw.stay_id.isin(active_stays)].copy()
    relevant.sort_values(
        ["stay_id","rel_hours"],kind="stable",inplace=True
    )
    result=stage2.loc[:,KEYS].copy().set_index("stay_id")
    grp=relevant.groupby("stay_id",sort=False)
    for source_col in ("vis_6drug_candidate","vis_7drug_candidate"):
        base=source_col.replace("_candidate","")
        result[base+"_pre6_max_candidate"]=grp[source_col].max().reindex(result.index)
        result[base+"_pre6_last_candidate"]=grp[source_col].last().reindex(result.index)
        result[base+"_pre6_mean_candidate"]=grp[source_col].mean().reindex(result.index)
    result["native_agent_present_anywhere"]=result.index.isin(active_stays).astype("int64")
    result["native_drug_bins_total"]=raw.groupby("stay_id").size().reindex(
        result.index,fill_value=0
    ).astype("int64")
    return result.reset_index(), {
        "native_drug_total_hourly_bins":len(raw),
        "native_drug_stays_with_positive_dose":len(active_stays),
        "native_drug_stays_with_any_hourly_bins":int(raw.stay_id.nunique()),
        "pre6_candidate_bins":len(relevant),
        "dose_source_units":{
            "vasopressin":"units/hour",
            "other_six_drugs":"micrograms/kg/min"
        },
        "candidate_score_formula_status":"NOT VALIDATED AGAINST STUDY-SPECIFIC ORIGINAL VIS"
    }
