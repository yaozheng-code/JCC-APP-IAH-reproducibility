# Methods and variables — local MIMIC-IV v3.1 study reproduction

## Study population and index time

Adult ICU stays (`admission_age >=18`) are identified from `mimiciv_derived.icustay_detail`. Intra-abdominal pressure uses `mimiciv_icu.chartevents` itemid 223763. The confirmation index time T0 is the first second consecutive IAP measurement >=12 mmHg, with measurements spaced no more than 12 h apart and in the valid (0,50] mmHg range. The original first adult hospital/ICU stay selection is applied by the cohort SQL.

## Exposure and outcomes

- `app_t0 = mbp_pre6_last - iap_t0` in mmHg.
- Primary binary exposure `APP<60 mmHg`; exploratory `APP<50 mmHg`.
- Continuous exposure: a 10 mmHg lower APP.
- 28-day mortality is based on death date relative to T0 calendar date (0–28 days).
- Hospital mortality is based on `hospital_expire_flag`.

## Main statistical adjustment fields

`admission_age`, `female`, `net_balance_pre24_L`, `albumin_pre24_last`, `vis_pre6_max`, `lactate_pre6_last`, `creat_pre24_last`, `k_pre24_last`, `sofa_dynamic_first`, `aki_dynamic_first`, `charlson`, `weight`.

The primary adjusted model additionally includes IAP at T0. The supplied unchanged `run_app_iah.py` handles missing covariates by within-analysis-sample median imputation and missingness indicator variables.

## Clinical definitions and time windows

- SOFA and AKI `dynamic_first` are the earliest eligible nonmissing ICU-period scores. First AKI **stage** is not synonymous with time of AKI onset; the latter is defined by first KDIGO stage >=1.
- Charlson uses the MIMIC derived admission-level comorbidity index.
- Weight is the earliest available record within the same hospital admission, recovered from the appropriate ICU stay in `weight_durations`.
- Sepsis and invasive ventilation follow the original SQL's per-stay record selection rules.
- T0-relative laboratory lookback windows and hourly fluid/urine bins use the definitions implemented in `native_stage1.py` and `native_time_window_transform.py`.
- Source-level eligibility conditions are applied to fluid input, output and urine to preserve genuinely unrecorded trajectories as missing instead of zero.

## Six-drug vasoactive-inotropic score

The investigator-specified workflow uses:

`VIS = dopamine + dobutamine + 100*epinephrine + 100*norepinephrine + 10*milrinone + 10000*vasopressin`

Vasopressin is converted from the MIMIC source unit of units/hour to units/kg/min by dividing by `60 * bodyweight_kg`; other drug terms are in µg/kg/min. Hourly scores use overlap-duration-weighted dose bins. The seven-drug exploratory score is not used as the manuscript model adjustment variable.

## Statistical outputs

The original unchanged `run_app_iah.py` generates main logistic models, tables by APP quartile, sensitivity and supplementary analyses, nonlinear adjusted curves, bootstrap discrimination metrics and 13 aggregate numerical benchmark checks.

`PASS` in the benchmark table means that an aggregate estimate met its specified tolerance; it should not be interpreted as a guarantee of exact patient-level equality. The package's aim is **independent approximate reproducibility**, with results and environmental details reported as obtained.

## Protection and dependencies

The software operates on the reproducer's local authorized MIMIC-IV v3.1 deployment. Do not distribute patient-level ICU/hospital IDs or raw extracted clinical data. The derived table versions and Python dependencies should be documented by each reproducer.
