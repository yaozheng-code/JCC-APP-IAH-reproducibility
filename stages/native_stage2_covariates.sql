-- JCC APP/IAH | MIMIC-IV v3.1 native clinical covariates, Stage 2.
-- READ ONLY. Input: array of 641 previously selected ICU stay IDs.
--
-- SCIENTIFIC DEFINITIONS (from the original study's SQL source records):
-- 1) sofa_dynamic_first: first non-NULL sofa_24hours by starttime in ICU stay.
-- 2) aki_dynamic_first: first non-NULL aki_stage by charttime in ICU stay.
-- 3) charlson: charlson_comorbidity_index for the same hadm_id.
-- 4) weight: first weight_durations record in the same HOSPITAL admission,
--    ordered by starttime (recovered via mimiciv_icu.icustays).
-- 5) sepsis: EXISTS any sepsis3 record for the stay (the original source
--    uses row existence, not a predicate on the sepsis3 boolean).
-- 6) InvasiveVent: EXISTS a ventilation_status='InvasiveVent' row anywhere
--    in the stay, with no extra restriction at T0.
-- Additional diagnostic sepsis3_true flag is NOT used in the manuscript.
-- Aggregated study scores sofa and aki_score mirror whole-stay maxima.
-- IMPORTANT: dynamic_first means first in ICU, NOT first after IAH T0.
--
-- Do not introduce a different time window, alter case definitions, or
-- compute VIS from arbitrary vasopressor-dose coefficients.

WITH c AS (
    SELECT d.subject_id, d.hadm_id, d.stay_id,
           d.icu_intime, d.icu_outtime
    FROM mimiciv_derived.icustay_detail d
    WHERE d.stay_id = ANY(%s::bigint[])
)
SELECT
    c.subject_id, c.hadm_id, c.stay_id,

    first_sofa.sofa_24hours AS sofa_dynamic_first,
    first_sofa.starttime AS sofa_dynamic_first_time,

    first_aki.aki_stage AS aki_dynamic_first,
    first_aki.charttime AS aki_dynamic_first_time,

    first_weight.weight AS weight,
    first_weight.starttime AS weight_first_time,
    first_weight.measurement_stay_id AS weight_source_stay_id,

    ch.charlson_comorbidity_index AS charlson,

    CASE WHEN EXISTS (
        SELECT 1 FROM mimiciv_derived.sepsis3 s
        WHERE s.stay_id = c.stay_id
    ) THEN 1 ELSE 0 END AS sepsis,

    CASE WHEN EXISTS (
        SELECT 1 FROM mimiciv_derived.sepsis3 s
        WHERE s.stay_id = c.stay_id AND s.sepsis3 IS TRUE
    ) THEN 1 ELSE 0 END AS sepsis3_true_diagnostic,

    CASE WHEN EXISTS (
        SELECT 1 FROM mimiciv_derived.ventilation v
        WHERE v.stay_id = c.stay_id
          AND v.ventilation_status = 'InvasiveVent'
    ) THEN 1 ELSE 0 END AS "InvasiveVent",

    score_max.sofa AS sofa,
    kdigo_max.aki_score AS aki_score

FROM c

LEFT JOIN LATERAL (
    SELECT s.starttime, s.sofa_24hours
    FROM mimiciv_derived.sofa s
    WHERE s.stay_id = c.stay_id
      AND s.sofa_24hours IS NOT NULL
      AND s.starttime >= c.icu_intime
      AND s.starttime <  c.icu_outtime
    ORDER BY s.starttime, s.endtime
    LIMIT 1
) first_sofa ON TRUE

LEFT JOIN LATERAL (
    SELECT k.charttime, k.aki_stage
    FROM mimiciv_derived.kdigo_stages k
    WHERE k.stay_id = c.stay_id
      AND k.aki_stage IS NOT NULL
      AND k.charttime >= c.icu_intime
      AND k.charttime <  c.icu_outtime
    ORDER BY k.charttime, k.aki_stage
    LIMIT 1
) first_aki ON TRUE

LEFT JOIN LATERAL (
    -- Original study 25270 SQL joined weight_durations by hadm_id;
    -- standard weight_durations has only stay_id, so use ICU stays to
    -- recover all weight records within the same hospital admission.
    SELECT w.weight, w.starttime, w.stay_id AS measurement_stay_id
    FROM mimiciv_icu.icustays i
    JOIN mimiciv_derived.weight_durations w
      ON w.stay_id = i.stay_id
    WHERE i.hadm_id = c.hadm_id
    ORDER BY w.starttime, w.endtime, w.stay_id
    LIMIT 1
) first_weight ON TRUE

LEFT JOIN LATERAL (
    SELECT x.charlson_comorbidity_index
    FROM mimiciv_derived.charlson x
    WHERE x.hadm_id = c.hadm_id
    ORDER BY x.subject_id
    LIMIT 1
) ch ON TRUE

LEFT JOIN LATERAL (
    SELECT MAX(s.sofa_24hours) AS sofa
    FROM mimiciv_derived.sofa s
    WHERE s.stay_id = c.stay_id
) score_max ON TRUE

LEFT JOIN LATERAL (
    SELECT MAX(k.aki_stage) AS aki_score
    FROM mimiciv_derived.kdigo_stages k
    WHERE k.stay_id = c.stay_id
) kdigo_max ON TRUE

ORDER BY c.stay_id;
