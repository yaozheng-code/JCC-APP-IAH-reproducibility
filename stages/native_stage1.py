#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JCC APP/IAH — LOCAL MIMIC-IV extraction, Stage 1 (native PostgreSQL only).

INPUT: 641-stay cohort generated from the user's own local MIMIC-IV database.
NO INPUT FROM THE 14 PREVIOUSLY EXPORTED FILES.
Uses 10 native read-only SQL source queries to reconstruct the IAP/T0 clock,
MAP/APP, laboratory/IO records, raw vasoactive drug dose timelines, and
death outcomes. Study VIS is intentionally NOT calculated until the original
weight normalization and formula are verified.

This is Stage 1, NOT the complete 12-covariate/13-benchmark analysis. Additional
weight, Charlson, invasive ventilation, sepsis, dynamic SOFA/AKI extraction and
verification will be done in Stage 2.

All patient-level output stays in the local authorized MIMIC environment.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from native_time_window_transform import (
    reconstruct_t0, matching_t0_trajectory, aggregate_window
)

HERE = Path(__file__).resolve().parent
INPUT_641 = (Path.home()/"mimic-project"/"results"/"jcc_app_iah"
             /"Final_cohort_641_LOCAL_ONLY.csv")
OUTPUT = (Path.home()/"mimic-project"/"results"/"jcc_app_iah_native_stage1")
DATABASE = "mimiciv"

# For installations with a different database name or connection string, edit
# DATABASE or pass --dsn via command line; NEVER include passwords in a paper ZIP.
TAGS = ["25252","25255","25256","25257","25258","25259",
        "25262","25263","25265","25267"]
FIELDS = {
    "25252":"custom_膀胱压力",
    "25255":"mbp", "25256":"lactate", "25257":"creatinine",
    "25258":"potassium", "25259":"albumin",
    "25262":"input_all", "25263":"output_all",
    "25265":"output_urine", "25267":"norepinephrine_mcgkgmin"
}
REQUIRED = {
    "mimiciv_derived.icustay_detail": [
        "subject_id","hadm_id","stay_id","icu_intime","icu_outtime",
        "admittime","dischtime","admission_age","gender","dod","hospital_expire_flag"
    ],
    "mimiciv_icu.chartevents":["hadm_id","charttime","itemid","valuenum"],
    "mimiciv_derived.vitalsign":["stay_id","charttime","mbp"],
    "mimiciv_derived.bg":["hadm_id","charttime","lactate","potassium","specimen"],
    "mimiciv_derived.chemistry":["hadm_id","charttime","creatinine","potassium","albumin"],
    "mimiciv_derived.urine_output":["stay_id","charttime","urineoutput"],
    "mimiciv_derived.vasoactive_agent":[
        "stay_id","starttime","endtime","dopamine","epinephrine",
        "norepinephrine","phenylephrine","vasopressin",
        "dobutamine","milrinone"
    ],
    "mimiciv_icu.inputevents":["hadm_id","starttime","endtime","amount","amountuom"],
    "mimiciv_icu.outputevents":["hadm_id","charttime","value"],
    "mimiciv_hosp.admissions":["subject_id","hadm_id"]
}
LATER_TABLES=[
    "mimiciv_derived.sofa","mimiciv_derived.kdigo_stages",
    "mimiciv_derived.charlson","mimiciv_derived.weight_durations",
    "mimiciv_derived.sepsis3","mimiciv_derived.ventilation",
]


def numeric_identifier_cohort(p:Path) -> pd.DataFrame:
    if not p.is_file():
        raise FileNotFoundError(
            f"找不到之前本地数据库生成的641例名单：{p}\n"
            "请先确认你保存的 Final_cohort_641_LOCAL_ONLY.csv。"
        )
    data=pd.read_csv(p,usecols=["subject_id","hadm_id","stay_id"])
    for c in ("subject_id","hadm_id","stay_id"):
        v=pd.to_numeric(data[c],errors="raise")
        if v.isna().any() or ((v%1)!=0).any():
            raise ValueError(f"{c}包含空值或无效编号")
        data[c]=v.astype("int64")
    if data.stay_id.duplicated().any() or len(data)!=641:
        raise ValueError(f"本地队列人数/唯一性异常: n={len(data)}")
    return data.sort_values("stay_id").reset_index(drop=True)


def check_schema(cursor):
    query="""
        SELECT a.attname
        FROM pg_attribute a
        WHERE a.attrelid = to_regclass(%s)
          AND a.attnum>0 AND NOT a.attisdropped
    """
    missing={}
    for table, columns in REQUIRED.items():
        cursor.execute(query,(table,))
        actual={x[0] for x in cursor.fetchall()}
        absent=sorted(set(columns)-actual)
        if absent:
            missing[table]=absent
    if missing:
        raise RuntimeError(
            "本地数据库与原研究提取SQL需要的表或字段不一致。\n"
            +json.dumps(missing,ensure_ascii=False,indent=2)
        )
    # Check additional concepts for next stage without requiring them here.
    status={}
    for table in LATER_TABLES:
        cursor.execute("SELECT to_regclass(%s) IS NOT NULL",(table,))
        status[table]=bool(cursor.fetchone()[0])
    print("第一阶段10类原始表/字段：检查通过",flush=True)
    print("第二阶段所需派生表已存在：",
          sum(status.values()),"/",len(status),flush=True)
    return status


def query_frame(cursor,sql:str,ids:list[int]):
    cursor.execute(sql,(ids,))
    columns=[x.name for x in cursor.description]
    rows=cursor.fetchall()
    return pd.DataFrame(rows,columns=columns)


def source_tables(cursor, ids:list[int]):
    sql="""
    SELECT d.subject_id,d.hadm_id,d.stay_id,d.admission_age,d.gender,
           d.dod,d.hospital_expire_flag,
           d.icu_intime,d.icu_outtime,d.admittime,d.dischtime
    FROM mimiciv_derived.icustay_detail d
    WHERE d.stay_id=ANY(%s::bigint[])
    ORDER BY d.stay_id
    """
    return query_frame(cursor,sql,ids)


def build_derived(first:pd.DataFrame, details:pd.DataFrame,
                  sources:dict[str,pd.DataFrame]) -> pd.DataFrame:
    df=details.merge(
        first[["stay_id","t0_time","t0_hr","iap_t0","iap_prev",
               "iap_gap_hr","iap_n_total","iap_n_ge12"]],
        on="stay_id",how="inner",validate="one_to_one"
    ).set_index("stay_id").sort_index()
    t0hr=df["t0_hr"]

    laboratory_rules=[
        ("25255",-6,0,"mbp_pre6",
            {"last":"last","last_rel":"last_rel","mean":"mean",
             "min":"min","max":"max","n":"n"}),
        ("25256",-6,0,"lactate_pre6",
            {"last":"last","last_rel":"last_rel","max":"max","n":"n"}),
        ("25257",-24,0,"creat_pre24",
            {"last":"last","last_rel":"last_rel","max":"max","n":"n"}),
        ("25257",0,24,"creat_0_24",
            {"last":"last","last_rel":"last_rel","max":"max","n":"n"}),
        ("25258",-24,0,"k_pre24",
            {"last":"last","last_rel":"last_rel","min":"min","n":"n"}),
        ("25258",0,24,"k_0_24",
            {"last":"last","last_rel":"last_rel","min":"min","n":"n"}),
        ("25259",-24,0,"albumin_pre24",
            {"last":"last","last_rel":"last_rel","n":"n"}),
        ("25259",-48,0,"albumin_pre48",
            {"last":"last","last_rel":"last_rel","n":"n"})
    ]
    rel={}
    for tag,table in sources.items():
        if tag=="25252":
            continue
        rel[tag]=matching_t0_trajectory(table,t0hr)

    for tag,start,end,prefix,mapping in laboratory_rules:
        a=aggregate_window(rel[tag],FIELDS[tag],start,end,inclusive="both")
        for suffix,stat in mapping.items():
            df[f"{prefix}_{suffix}"]=(a[stat].reindex(df.index)
                                     if not a.empty else np.nan)
    # The source study provided VIS directly. Local standard vasoactive_agent
    # contains drug doses, NOT VIS, and vasopressin has units/hour. Never
    # replace manuscript VIS with an unverified formula or zero.
    for suffix in ("last","last_rel","mean","max","n"):
        df[f"vis_pre6_{suffix}"] = np.nan

    for prefix,tag in [("input","25262"),("output","25263"),
                       ("urine","25265")]:
        for start,end,label in [(-24,0,"pre24"),(0,24,"0_24"),
                                 (24,48,"24_48")]:
            a=aggregate_window(rel[tag],FIELDS[tag],start,end,inclusive="left")
            df[f"{prefix}_{label}"]=(a["sum"].reindex(df.index)
                                    if not a.empty else np.nan)

    for label in ["pre24","0_24","24_48"]:
        df[f"net_balance_{label}"]=df[f"input_{label}"]-df[f"output_{label}"]
    df["net_balance_pre24_L"]=df["net_balance_pre24"]/1000
    df["female"]=(df["gender"].astype(str).str.upper()=="F").astype(int)

    df["app_t0"]=df["mbp_pre6_last"]-df["iap_t0"]
    df["app_low60"]=np.where(df.app_t0.notna(),
                              (df.app_t0<60).astype(int),np.nan)
    df["app_low50"]=np.where(df.app_t0.notna(),
                              (df.app_t0<50).astype(int),np.nan)
    df["app_min_pre6"]=df["mbp_pre6_min"]-df["iap_t0"]
    df["lower_app_per10"]=-df["app_t0"]/10
    df["lower_app_min_per10"]=-df["app_min_pre6"]/10
    df["app_min_low60"]=np.where(df.app_min_pre6.notna(),
                                  (df.app_min_pre6<60).astype(int),np.nan)
    df["sepsis"]=np.nan  # Reserved for independently verified Stage 2 SQL.
    df["InvasiveVent"]=np.nan
    df["sofa_dynamic_first"]=np.nan
    df["aki_dynamic_first"]=np.nan
    df["charlson"]=np.nan
    df["weight"]=np.nan
    # The date of death has calendar-date resolution, as in the original pipeline.
    dod=pd.to_datetime(df.dod,errors="coerce").dt.normalize()
    t0=pd.to_datetime(df.t0_time,errors="coerce").dt.normalize()
    df["death_28d_from_t0"]=(dod-t0).dt.days.between(0,28).astype(int)
    df["hospital_expire_flag"]=pd.to_numeric(
        df.hospital_expire_flag,errors="coerce"
    )
    df["t0_to_icuout_hr"]=(
        pd.to_datetime(df.icu_outtime,errors="coerce")
        -pd.to_datetime(df.t0_time,errors="coerce")
    ).dt.total_seconds()/3600
    if df.app_t0.notna().all():
        df["app_q"]=pd.qcut(
            df.app_t0,q=4,
            labels=["Q1 lowest","Q2","Q3","Q4 highest"]
        )
    # Make clear that this is NOT yet the completed analysis dataset.
    return df.reset_index()


def main():
    import argparse
    try:
        import psycopg
    except ImportError as exc:
        raise RuntimeError(
            "缺少psycopg，先运行 %pip install 'psycopg[binary]'"
        ) from exc

    p=argparse.ArgumentParser()
    p.add_argument("--database",default=DATABASE)
    p.add_argument("--dsn",default=None)
    p.add_argument("--cohort",type=Path,default=INPUT_641)
    p.add_argument("--outdir",type=Path,default=OUTPUT)
    args=p.parse_args()

    cohort=numeric_identifier_cohort(args.cohort.expanduser())
    ids=cohort.stay_id.astype("int64").tolist()
    output=args.outdir.expanduser()
    output.mkdir(parents=True,exist_ok=True)
    extracts=output/"native_source_trajectories"
    extracts.mkdir(parents=True,exist_ok=True)

    print("="*65,flush=True)
    print("JCC APP–IAH | 原始本地MIMIC数据库直接提取：第一阶段",flush=True)
    print("输入：先前本地PostgreSQL筛出的641例；不读取14份历史导出",flush=True)
    print("="*65,flush=True)

    with psycopg.connect(
        args.dsn if args.dsn else f"dbname={args.database}"
    ) as connection:
        with connection.cursor() as cursor:
            cursor.execute("SET TRANSACTION READ ONLY")
            cursor.execute("SELECT current_database()")
            dbname=cursor.fetchone()[0]
            print("数据库连接：",dbname,flush=True)
            print("PostgreSQL NUMERIC/Decimal → float64 时间窗兼容处理：已启用",flush=True)
            future_tables=check_schema(cursor)
            print("正在预检查10份SQL（EXPLAIN，无真实执行）...",flush=True)
            for preflight_tag in TAGS:
                preflight_sql=(HERE/f"native_extract_{preflight_tag}.sql").read_text(
                    encoding="utf-8")
                try:
                    cursor.execute("EXPLAIN " + preflight_sql, (ids,))
                except Exception as ex:
                    raise RuntimeError(
                        f"提取SQL {preflight_tag} 尚不兼容本地表结构：{ex}"
                    ) from ex
            print("10/10 SQL结构预检查通过",flush=True)
            details=source_tables(cursor,ids)
            if len(details)!=641 or details.stay_id.duplicated().any():
                raise RuntimeError("icustay_detail与641名单无法一对一对应")
            if set(details.stay_id)!=set(cohort.stay_id):
                raise RuntimeError("本地队列ID与底层ICU明细不一致")

            sources={}
            stats={}
            for tag in TAGS:
                sql=(HERE/f"native_extract_{tag}.sql").read_text(encoding="utf-8")
                print(f"直接查询MIMIC数据库：{tag} / {FIELDS[tag]} ...",
                      flush=True)
                tab=query_frame(cursor,sql,ids)
                if not {"stay_id","t_hours_from_icu_intime",FIELDS[tag]}.issubset(tab):
                    raise RuntimeError(f"{tag}查询缺少必需列：{tab.columns.tolist()}")
                for key in ("subject_id","hadm_id","stay_id"):
                    tab[key]=pd.to_numeric(tab[key],errors="raise").astype("int64")
                sources[tag]=tab
                stats[tag]=len(tab)
                print(f"    提取记录数：{len(tab)}",flush=True)

    # Reuse the validated time-window and T0 reconstruction implementations.
    first,_=reconstruct_t0(sources["25252"])
    print("本地原始IAP确认的目标ICU住院：",len(first),"/641",flush=True)
    if first.stay_id.duplicated().any() or (
        set(first.stay_id)!=set(cohort.stay_id)
    ):
        raise RuntimeError(
            "使用本地原始IAP重建T0时，确认的641名单与已选队列不一致。"
            "停止执行；不能使用历史CSV来填补或代替本地结果。"
        )

    result=build_derived(first,details,sources)
    if len(result)!=641 or result.stay_id.duplicated().any():
        raise RuntimeError("重建后不是641条唯一ICU住院记录")

    # Save only after the 641 T0 consistency checks have passed.
    for tag,tab in sources.items():
        tab.to_csv(
            extracts/f"native_{tag}_{FIELDS[tag]}.csv.gz",
            compression="gzip",index=False
        )
    outfile=output/"native_stage1_641_clinical_partial_LOCAL_ONLY.csv"
    result.to_csv(outfile,index=False)

    missing_app=int(result.app_t0.isna().sum())
    obs60=(int(result.app_low60.sum())
           if missing_app==0 else None)
    obs50=(int(result.app_low50.sum())
           if missing_app==0 else None)

    report={
        "run_time_utc":datetime.now(timezone.utc).isoformat(),
        "database":dbname,
        "sole_clinical_source":"Local authorized MIMIC-IV PostgreSQL v3.1",
        "prior_14_export_files_used_as_input":False,
        "cohort_stays":len(cohort),
        "T0_confirmed_stays":len(first),
        "source_event_rows":stats,
        "missing_APP_at_T0":missing_app,
        "APP_below_60":obs60,
        "APP_below_50":obs50,
        "28_day_deaths_from_T0":int(result.death_28d_from_t0.sum()),
        "hospital_deaths":int(result.hospital_expire_flag.sum()),
        "later_derived_table_presence":future_tables,
        "vasoactive_agent_dose_source":"mimiciv_derived.vasoactive_agent",
        "vis_score_status":"NOT COMPUTED: study-specific VIS formula and vasopressin weight conversion not yet verified",
        "not_yet_completed_columns":[
            "sofa_dynamic_first","aki_dynamic_first","charlson",
            "weight","sepsis","InvasiveVent","vis_pre6_last",
            "vis_pre6_mean","vis_pre6_max"
        ],
        "manuscript_models_recomputed":False,
        "scope":"Stage 1 of native-source clinical reconstruction; no model claims.",
    }
    (output/"native_stage1_summary.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8"
    )

    print("\n"+"="*65,flush=True)
    print("本地数据库临床数据提取结果（第一阶段）",flush=True)
    print("="*65,flush=True)
    print("唯一ICU住院：",len(result),flush=True)
    print("本地T0确认：",len(first),flush=True)
    print("APP缺失：",missing_app,flush=True)
    print("APP<60：",obs60,"（历史参考367；不强制修改）",flush=True)
    print("APP<50：",obs50,"（历史参考180；不强制修改）",flush=True)
    print("T0后28天死亡：",report["28_day_deaths_from_T0"],
          "（历史参考225）",flush=True)
    print("院内死亡：",report["hospital_deaths"],"（历史参考225）",flush=True)
    print("尿量：已按本地stay_id结构重新提取",flush=True)
    print("血管活性药：已提取7种药物的每小时剂量（VIS总分暂不推算）",flush=True)
    print("已保存本地文件：",output,flush=True)
    print("下一阶段：核实VIS及补齐首次SOFA/AKI、Charlson、体重、脓毒症、机械通气。",flush=True)
    print("这一阶段不运行统计模型，也不读取历史14份导出。",flush=True)
    print("="*65,flush=True)


if __name__=="__main__":
    main()
