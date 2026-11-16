-- MIMIC-IV v3.1 PostgreSQL, native vasoactive-agent source extraction.
-- The original study's source table supplied precomputed VIS (vis_score).
-- The standard MIT-LCP mimiciv_derived.vasoactive_agent table DOES NOT supply VIS;
-- its vasopressin rate is units/hour, other agents are mcg/kg/min.
-- To preserve model fidelity, output the seven raw time-weighted drug doses and
-- defer VIS recomputation until the original formula and weight handling are verified.
-- This is NOT a substitute for the original VIS source.
-- READ ONLY; 641 authorized ICU stay IDs are passed as a single parameter.
WITH traj_cohort AS (
    SELECT d.subject_id, d.hadm_id, d.stay_id, d.icu_intime, d.icu_outtime
    FROM mimiciv_derived.icustay_detail d
    WHERE d.stay_id = ANY(%s::bigint[])
),
bin_grid AS (
    SELECT c.subject_id, c.hadm_id, c.stay_id, c.icu_intime,
           b.bin_idx,
           c.icu_intime + b.bin_idx * INTERVAL '1 hour' AS bin_start,
           LEAST(c.icu_intime + (b.bin_idx+1)*INTERVAL '1 hour',
                 c.icu_outtime) AS bin_end
    FROM traj_cohort c
    CROSS JOIN LATERAL generate_series(
        0,
        GREATEST(0, CEIL(EXTRACT(EPOCH FROM (c.icu_outtime-c.icu_intime))
                         / 3600.0)::int - 1)
    ) b(bin_idx)
),
hourly_doses AS (
    SELECT g.stay_id, g.bin_idx,
           SUM(COALESCE(p.dopamine,0) * EXTRACT(EPOCH FROM
               (LEAST(p.endtime,g.bin_end)-GREATEST(p.starttime,g.bin_start))))
             / NULLIF(MAX(EXTRACT(EPOCH FROM(g.bin_end-g.bin_start))),0)
             AS dopamine_mcgkgmin,
           SUM(COALESCE(p.epinephrine,0) * EXTRACT(EPOCH FROM
               (LEAST(p.endtime,g.bin_end)-GREATEST(p.starttime,g.bin_start))))
             / NULLIF(MAX(EXTRACT(EPOCH FROM(g.bin_end-g.bin_start))),0)
             AS epinephrine_mcgkgmin,
           SUM(COALESCE(p.norepinephrine,0) * EXTRACT(EPOCH FROM
               (LEAST(p.endtime,g.bin_end)-GREATEST(p.starttime,g.bin_start))))
             / NULLIF(MAX(EXTRACT(EPOCH FROM(g.bin_end-g.bin_start))),0)
             AS norepinephrine_mcgkgmin,
           SUM(COALESCE(p.phenylephrine,0) * EXTRACT(EPOCH FROM
               (LEAST(p.endtime,g.bin_end)-GREATEST(p.starttime,g.bin_start))))
             / NULLIF(MAX(EXTRACT(EPOCH FROM(g.bin_end-g.bin_start))),0)
             AS phenylephrine_mcgkgmin,
           SUM(COALESCE(p.vasopressin,0) * EXTRACT(EPOCH FROM
               (LEAST(p.endtime,g.bin_end)-GREATEST(p.starttime,g.bin_start))))
             / NULLIF(MAX(EXTRACT(EPOCH FROM(g.bin_end-g.bin_start))),0)
             AS vasopressin_units_hour,
           SUM(COALESCE(p.dobutamine,0) * EXTRACT(EPOCH FROM
               (LEAST(p.endtime,g.bin_end)-GREATEST(p.starttime,g.bin_start))))
             / NULLIF(MAX(EXTRACT(EPOCH FROM(g.bin_end-g.bin_start))),0)
             AS dobutamine_mcgkgmin,
           SUM(COALESCE(p.milrinone,0) * EXTRACT(EPOCH FROM
               (LEAST(p.endtime,g.bin_end)-GREATEST(p.starttime,g.bin_start))))
             / NULLIF(MAX(EXTRACT(EPOCH FROM(g.bin_end-g.bin_start))),0)
             AS milrinone_mcgkgmin
    FROM bin_grid g
    JOIN mimiciv_derived.vasoactive_agent p
      ON p.stay_id=g.stay_id
     AND p.starttime < g.bin_end
     AND p.endtime > g.bin_start
     AND p.endtime > p.starttime
    GROUP BY g.stay_id, g.bin_idx
)
SELECT g.subject_id, g.hadm_id, g.stay_id,
       ROUND(g.bin_idx::numeric,2) AS "t_hours_from_icu_intime",
       ROUND(COALESCE(h.dopamine_mcgkgmin,0)::numeric,6) AS dopamine_mcgkgmin,
       ROUND(COALESCE(h.epinephrine_mcgkgmin,0)::numeric,6) AS epinephrine_mcgkgmin,
       ROUND(COALESCE(h.norepinephrine_mcgkgmin,0)::numeric,6) AS norepinephrine_mcgkgmin,
       ROUND(COALESCE(h.phenylephrine_mcgkgmin,0)::numeric,6) AS phenylephrine_mcgkgmin,
       ROUND(COALESCE(h.vasopressin_units_hour,0)::numeric,6) AS vasopressin_units_hour,
       ROUND(COALESCE(h.dobutamine_mcgkgmin,0)::numeric,6) AS dobutamine_mcgkgmin,
       ROUND(COALESCE(h.milrinone_mcgkgmin,0)::numeric,6) AS milrinone_mcgkgmin
FROM bin_grid g
LEFT JOIN hourly_doses h USING(stay_id,bin_idx)
ORDER BY g.stay_id,g.bin_idx;
