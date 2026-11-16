# Independent reproducer's local run report (fill in locally)

**No patient-level CSVs or stay IDs should be attached to this report.**

- MIMIC-IV database version: 3.1 / ______
- PostgreSQL version: ______
- Derived concept scripts version/commit/date: ______
- Python version: ______
- Code ZIP SHA256: ______
- Run UTC date: ______
- Number of IAH candidate ICU stays: ______ (reference 819)
- Final eligible first adult ICU stays: ______ (reference 641)
- T0 confirmations from native IAP: ______ (reference 641)
- APP<60: ______ (reference 367)
- APP<50: ______ (reference 180)
- T0+28-day deaths: ______ (reference 225)
- Six-drug VIS observed / missing: ______ / ______ (prior local 469 / 172)
- Source-ineligible fluid input / urine: ______ / ______ (prior local 26 / 4)
- Adjusted APP<60 vs 28-day mortality OR: ______ (published reference approximately 2.03)
- OR 95% CI, p-value: ______
- Aggregate benchmark comparisons (PASS/REVIEW): ______ / ______
- Source-level discrepancies or errors: ______
- Files checked for inadvertent patient data release: Yes / No
- Reproducer conclusion: independent approximate / review needed / not reproduced
- Notes: "Aggregate PASS does not mean individual records are identical."

Do not alter a patient-level variable to reproduce a published OR. Report differences explicitly.
