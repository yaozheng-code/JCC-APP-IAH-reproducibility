-- Original IAH-positive candidate cohort, from previously tested local MIMIC-IV 3.1 query.
-- Read-only; no historical source exports.
WITH
adult_icu AS MATERIALIZED (
    SELECT subject_id, hadm_id, stay_id,
           admission_age, admittime, icu_intime, icu_outtime
    FROM mimiciv_derived.icustay_detail
    WHERE admission_age >= 18
      AND icu_intime IS NOT NULL
      AND icu_outtime IS NOT NULL
),
stays_with_iap AS (
    SELECT s.*
    FROM adult_icu s
    WHERE EXISTS (
        SELECT 1
        FROM mimiciv_icu.chartevents ce
        WHERE ce.hadm_id = s.hadm_id
          AND ce.itemid = 223763
          AND ce.valuenum IS NOT NULL
          AND ce.charttime >= s.icu_intime
          AND ce.charttime < s.icu_outtime
    )
),
ranked AS (
    SELECT s.*,
           DENSE_RANK() OVER (
               PARTITION BY subject_id ORDER BY admittime
           ) AS filtered_hospital_seq,
           DENSE_RANK() OVER (
               PARTITION BY hadm_id ORDER BY icu_intime
           ) AS filtered_icu_seq
    FROM stays_with_iap s
),
initial_selection AS (
    SELECT *
    FROM ranked
    WHERE filtered_hospital_seq = 1
      AND filtered_icu_seq = 1
),
iap_valid AS (
    SELECT s.subject_id, s.hadm_id, s.stay_id,
           ce.charttime AS measurement_time,
           ce.valuenum AS iap
    FROM initial_selection s
    JOIN mimiciv_icu.chartevents ce
      ON ce.hadm_id = s.hadm_id
     AND ce.itemid = 223763
     AND ce.valuenum > 0
     AND ce.valuenum <= 50
     AND ce.charttime >= s.icu_intime
     AND ce.charttime < s.icu_outtime
),
iap_in_order AS (
    SELECT e.*,
           LAG(iap) OVER (
               PARTITION BY stay_id ORDER BY measurement_time
           ) AS prior_iap,
           LAG(measurement_time) OVER (
               PARTITION BY stay_id ORDER BY measurement_time
           ) AS prior_measurement_time
    FROM iap_valid e
),
iah_confirmed AS (
    SELECT DISTINCT subject_id, hadm_id, stay_id
    FROM iap_in_order
    WHERE iap >= 12
      AND prior_iap >= 12
      AND measurement_time >= prior_measurement_time
      AND measurement_time - prior_measurement_time <= INTERVAL '12 hours'
)
SELECT subject_id, hadm_id, stay_id
FROM iah_confirmed
ORDER BY stay_id;
