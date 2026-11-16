-- MIMIC-IV v3.1 | local PostgreSQL extraction | patient-level, credentialed use ONLY
-- MIMIC-IV native urine_output is stay_id-keyed, not hadm_id-keyed.
-- Scope limited to the already confirmed local 641 ICU stay_id set.
-- READ ONLY; no custom timestamp functions or preview LIMITs.
WITH traj_cohort AS (
 SELECT d.subject_id, d.hadm_id, d.stay_id,
        d.icu_intime, d.icu_outtime, d.admittime, d.dischtime
 FROM mimiciv_derived.icustay_detail d
 WHERE d.stay_id = ANY(%s::bigint[])
)
, bin_grid as (
    select c.subject_id, c.hadm_id, c.stay_id,
           c.icu_intime as win_start,
           c.icu_outtime   as win_end,
           b.bin_idx,
           c.icu_intime
               + b.bin_idx * interval '1 hours' as bin_start,
           least(
               c.icu_intime
                   + (b.bin_idx + 1) * interval '1 hours',
               c.icu_outtime
           ) as bin_end
    from traj_cohort c
    cross join lateral generate_series(
        0,
        greatest(0, ceil(extract(epoch from (
            c.icu_outtime
            - c.icu_intime
        )) / (1 * 3600.0))::int - 1)
    ) as b(bin_idx)
)
, overlap as (
    select g.stay_id, g.bin_idx, sum(p.urineoutput) as value_in_bin
    from bin_grid g
    join mimiciv_derived.urine_output p
      on p.stay_id = g.stay_id
     and p.charttime >= g.bin_start
     and p.charttime <  g.bin_end
    where p.urineoutput is not null
      
    group by g.stay_id, g.bin_idx
)
select g.subject_id,
       g.hadm_id,
       g.stay_id,
       round(
           (g.bin_idx * extract(epoch from interval '1 hours')
            / extract(epoch from interval '1 hour'))::numeric,
           2
       ) as "t_hours_from_icu_intime",
       round(coalesce(o.value_in_bin, 0)::numeric, 1) as "output_urine"
from bin_grid g
left join overlap o using (stay_id, bin_idx)
order by g.stay_id, g.bin_idx
;
