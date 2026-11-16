#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Original study T0 and relative-window transformation with a type fix.

The IAP/T0 definitions, event ordering, and time-window boundaries are unchanged.
Only PostgreSQL numeric/Decimal values are explicitly cast to float64 before
arithmetic so psycopg3 and CSV inputs execute consistently.
"""
from __future__ import annotations
import pandas as pd
import numpy as np

KEYS = ["subject_id","hadm_id","stay_id"]

def reconstruct_t0(iap_raw: pd.DataFrame):
    source=iap_raw.copy()
    source["event_time"]=pd.to_datetime(source["event_time"])
    # PostgreSQL round(numeric, n) is returned as decimal.Decimal by psycopg3.
    # Normalize all time and pressure numerics BEFORE doing arithmetic so that
    # both CSV fixtures and native PostgreSQL use identical numeric types.
    source["t_hours_from_icu_intime"]=pd.to_numeric(
        source["t_hours_from_icu_intime"],errors="coerce"
    ).astype("float64")
    source["custom_膀胱压力"]=pd.to_numeric(
        source["custom_膀胱压力"],errors="coerce"
    ).astype("float64")
    valid=source[source["custom_膀胱压力"].between(0,50,inclusive="right")].copy()
    valid=valid.sort_values(["stay_id","event_time"],kind="stable")
    valid["iap_prev"]=valid.groupby("stay_id")["custom_膀胱压力"].shift()
    valid["time_prev"]=valid.groupby("stay_id")["event_time"].shift()
    valid["prev_hour"]=valid.groupby("stay_id")["t_hours_from_icu_intime"].shift()
    valid["gap_from_event_time"]=(valid.event_time-valid.time_prev).dt.total_seconds()/3600
    qualifying=(valid["custom_膀胱压力"]>=12)&(valid.iap_prev>=12)&valid.gap_from_event_time.between(0,12)
    first=valid[qualifying].drop_duplicates("stay_id",keep="first").copy()
    first.rename(columns={
        "event_time":"t0_time",
        "t_hours_from_icu_intime":"t0_hr",
        "custom_膀胱压力":"iap_t0"
    },inplace=True)
    first["iap_gap_hr"]=first["t0_hr"]-first["prev_hour"]
    counts=valid.groupby("stay_id")["custom_膀胱压力"].agg(
        iap_n_total="size",
        iap_n_ge12=lambda x:int(x.ge(12).sum())
    )
    first=first.join(counts,on="stay_id")
    fields=KEYS+["t0_time","t0_hr","iap_t0","iap_prev","iap_gap_hr",
                 "iap_n_total","iap_n_ge12"]
    return first[fields].copy(),valid


def matching_t0_trajectory(frame, t0):
    if frame.empty:
        return frame.copy()
    result=frame.copy()
    result["t_hours_from_icu_intime"]=pd.to_numeric(
        result["t_hours_from_icu_intime"],errors="coerce"
    ).astype("float64")
    # t0_hr may originate from a PostgreSQL numeric/decimal.Decimal, too.
    result["t0_hr"]=pd.to_numeric(
        result["stay_id"].map(t0),errors="coerce"
    ).astype("float64")
    result["rel_hours"]=result["t_hours_from_icu_intime"]-result["t0_hr"]
    return result


def aggregate_window(raw:pd.DataFrame, field:str, start:float,end:float,
                     inclusive="both")->pd.DataFrame:
    if raw.empty:
        return pd.DataFrame()
    mask=raw.rel_hours.between(start,end,inclusive=inclusive)
    a=raw.loc[mask,["stay_id","rel_hours",field]].copy()
    if a.empty:
        return pd.DataFrame()
    a[field]=pd.to_numeric(a[field],errors="coerce").astype("float64")
    a.sort_values(["stay_id","rel_hours"],kind="stable",inplace=True)
    g=a.groupby("stay_id",sort=False)
    table=pd.DataFrame(index=g.size().index)
    table["last"]=g[field].last()
    table["last_rel"]=g["rel_hours"].last()
    table["max"]=g[field].max()
    table["min"]=g[field].min()
    table["mean"]=g[field].mean()
    table["n"]=g[field].count()
    table["sum"]=g[field].sum()
    return table

