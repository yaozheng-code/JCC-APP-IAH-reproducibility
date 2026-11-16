-- Select first hospital and ICU encounter among adult MIMIC-IV stays.
-- One bound array of candidate stay_id values; no PHI bundled.
WITH adult_icu AS (
    SELECT
        subject_id, hadm_id, stay_id,
        admittime, icu_intime,
        DENSE_RANK() OVER (
            PARTITION BY subject_id ORDER BY admittime
        ) AS first_filtered_hospital,
        DENSE_RANK() OVER (
            PARTITION BY hadm_id ORDER BY icu_intime
        ) AS first_filtered_icu
    FROM mimiciv_derived.icustay_detail
    WHERE admission_age >= 18
),
original_first_stays AS (
    SELECT subject_id, hadm_id, stay_id
    FROM adult_icu
    WHERE first_filtered_hospital = 1
      AND first_filtered_icu = 1
      AND icu_intime IS NOT NULL
)
SELECT subject_id, hadm_id, stay_id
FROM original_first_stays
WHERE stay_id = ANY(%s::bigint[])
ORDER BY stay_id;
