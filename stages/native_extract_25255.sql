-- MIMIC-IV v3.1 | local PostgreSQL extraction | patient-level, credentialed use ONLY
-- Derived from the ORIGINAL study source-extraction query's event / hourly-bin logic.
-- Scope limited to the already confirmed local 641 ICU stay_id set.
-- READ ONLY; no custom timestamp functions or preview LIMITs.
WITH traj_cohort AS (
 SELECT d.subject_id, d.hadm_id, d.stay_id,
        d.icu_intime, d.icu_outtime, d.admittime, d.dischtime
 FROM mimiciv_derived.icustay_detail d
 WHERE d.stay_id = ANY(%s::bigint[])
)
, traj_events as (
        select
            mimiciv_derived.vitalsign.stay_id as join_key,
            mimiciv_derived.vitalsign.charttime as event_time,
            mimiciv_derived.vitalsign.mbp as value_raw
        from mimiciv_derived.vitalsign
        where mimiciv_derived.vitalsign.mbp is not null
            
)

select
    traj_cohort.subject_id,
    traj_cohort.hadm_id,
    traj_cohort.stay_id,
    round(
        (extract(epoch from (traj_events.event_time - traj_cohort.icu_intime))
         / extract(epoch from interval '1 hour'))::numeric,
        2
    ) as "t_hours_from_icu_intime",
    traj_events.event_time,
    traj_events.value_raw as "mbp"
from traj_cohort
join traj_events
    on traj_events.join_key = traj_cohort.stay_id
where traj_events.event_time >= traj_cohort.icu_intime
  and traj_events.event_time <  traj_cohort.icu_outtime
order by traj_cohort.stay_id, traj_events.event_time
;
