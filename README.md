# Local MIMIC-IV v3.1 reproducibility code — APP at confirmed IAH

This repository provides a **code-only independent approximate reproduction** workflow for the study *Abdominal perfusion pressure at confirmation of intra-abdominal hypertension and subsequent mortality in critically ill patients: a retrospective cohort study*.

Data come exclusively from the reproducer's authorized **local MIMIC-IV v3.1 PostgreSQL installation**. No patient-level data, historical exports, preselected cohorts, or original analysis CSV files are included or required.

## Scope

The pipeline independently identifies the adult IAH cohort and confirmation index time T0, reconstructs APP and laboratory/physiological variables, obtains clinical covariates from standard MIMIC derived tables, calculates the investigator-specified six-drug vasoactive-inotropic score, preserves measurement-source missingness, and runs the original unchanged manuscript statistical implementation.

Outputs include core models, sensitivity analyses, descriptive tables, figures, bootstrap AUCs, and thirteen aggregate benchmark comparisons.

## Prerequisites

- Credentialed access to MIMIC-IV **v3.1** in PostgreSQL.
- Schemas `mimiciv_icu`, `mimiciv_hosp`, and `mimiciv_derived`.
- Required derived concepts already installed: `icustay_detail`, `vitalsign`, `bg`, `chemistry`, `urine_output`, `vasoactive_agent`, `sofa`, `kdigo_stages`, `charlson`, `weight_durations`, `sepsis3`, and `ventilation`.
- Python 3.11+, a PostgreSQL user with SELECT permissions, and the dependencies in `requirements.txt`.

No credentials are stored in the package. Configure libpq environment variables or a local `.pgpass` file if needed.

## Run

In the extracted repository root:

```bash
python -m pip install -r requirements.txt
python run_all.py --database mimiciv --bootstrap 500 --seed 42
```

For a faster preliminary run, `--bootstrap 200` can be used. Specify a local private results directory with `--outdir`.

For Jupyter, see `Run_One_Click_Jupyter.ipynb` and `JUPYTER_一格启动_示例.txt`.

## Outputs and interpretation

The default output directory is:

`~/mimic-project/results/jcc_app_iah_native_approx_reproduction/`

It contains local clinical extracts and model inputs, manuscript tables, figures, and `REPRODUCTION_SUMMARY_AGGREGATE_ONLY.json`.

**Keep all patient-level records, stay identifiers, and extracted time series within your authorized MIMIC environment.** Review aggregate outputs before sharing, in accordance with your data use agreement.

Reproduction should be interpreted at the levels of study definitions, cohort construction, model specification, major estimates and uncertainty. Minor numerical variation across derived concept versions and local environments is possible. The benchmark audit reports observed results rather than forcing target values.

The investigator has successfully completed an end-to-end run of the integrated pipeline against a locally deployed MIMIC-IV v3.1 database. This confirms execution in the investigator's local environment; it does not guarantee that the workflow will run without adjustment on every computer or database installation. Successful reproduction depends on compatible database and derived-concept schemas, Python dependencies, PostgreSQL configuration, access permissions, and available computing resources. Users should document their environments, review execution logs, and verify the generated results.

Please cite the MIMIC-IV v3.1 dataset and applicable derived-concept sources according to PhysioNet requirements.
