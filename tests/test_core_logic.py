#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data-free validation: 641 synthetic ICU stays only.

No database connection, actual IDs, or patient records are involved.
Run from package root:
    python -m unittest discover -s tests -v
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

STAGES=Path(__file__).resolve().parents[1]/"stages"
sys.path.insert(0,str(STAGES))
from compute_vis import candidate_hourly_drug_scores, DRUGS
from source_eligibility import restore_source_missingness


class NativeStudySourceTests(unittest.TestCase):
    def setUp(self):
        n=641
        ids=pd.DataFrame({
            "stay_id":np.arange(100001,100001+n,dtype=np.int64),
            "hadm_id":np.arange(200001,200001+n,dtype=np.int64),
            "subject_id":np.arange(300001,300001+n,dtype=np.int64),
        })
        self.ids=ids
        self.native=ids.copy()
        self.native["t0_hr"]=10.0
        self.native["weight"]=75.0
        self.native["app_t0"]=60.0
        self.native["iap_t0"]=12.0
        self.native["mbp_pre6_last"]=72.0
        self.native["death_28d_from_t0"]=0
        self.native["hospital_expire_flag"]=0
        self.native["aki_dynamic_first"]=0
        self.native["sofa_dynamic_first"]=4
        self.native["vis_pre6_max"]=np.nan
        for label in ("pre24","0_24","24_48"):
            for source in ("input","output","urine"):
                self.native[f"{source}_{label}"]=100.0
            self.native[f"net_balance_{label}"]=0.0
        self.native["net_balance_pre24_L"]=0.0

    def test_VIS_six_drug_units_and_not_seven_drug(self):
        bins=self.ids.copy()
        bins["t_hours_from_icu_intime"]=9.0
        for name in DRUGS:
            bins[name]=0.0
        mask=bins.stay_id.eq(int(self.ids.stay_id.iloc[0]))
        bins.loc[mask,"dopamine_mcgkgmin"]=2.
        bins.loc[mask,"dobutamine_mcgkgmin"]=4.
        bins.loc[mask,"epinephrine_mcgkgmin"]=0.04
        bins.loc[mask,"norepinephrine_mcgkgmin"]=0.08
        bins.loc[mask,"milrinone_mcgkgmin"]=0.20
        bins.loc[mask,"phenylephrine_mcgkgmin"]=0.30
        bins.loc[mask,"vasopressin_units_hour"]=2.25
        scores,_=candidate_hourly_drug_scores(self.native,bins)
        # At 75 kg, vasopressin contributes 10000*(2.25/(75*60)) = 5.
        six=2+4+100*0.04+100*0.08+10*0.20+10000*2.25/(75*60)
        seven=six+10*0.30
        row=scores.loc[mask].iloc[0]
        self.assertAlmostEqual(row.vis_6drug_pre6_max_candidate,25.0,places=9)
        self.assertAlmostEqual(row.vis_6drug_pre6_max_candidate,six,places=9)
        self.assertAlmostEqual(row.vis_7drug_pre6_max_candidate,seven,places=9)
        self.assertEqual(int(scores.vis_6drug_pre6_max_candidate.notna().sum()),1)
        self.assertEqual(int(scores.vis_6drug_pre6_max_candidate.isna().sum()),640)

    def test_original_source_eligibility_restores_true_missing(self):
        flags=self.ids.copy()
        flags["has_input_source_25262"]=True
        flags["has_output_source_25263"]=True
        flags["has_urine_source_25265"]=True
        flags.loc[:25,"has_input_source_25262"]=False
        flags.loc[:3,"has_urine_source_25265"]=False
        result,summary=restore_source_missingness(self.native,flags)
        self.assertEqual(summary["input"]["source_ineligible_stays"],26)
        self.assertEqual(summary["urine"]["source_ineligible_stays"],4)
        self.assertEqual(int(result.net_balance_pre24_L.isna().sum()),26)
        self.assertEqual(int(result.urine_0_24.isna().sum()),4)
        self.assertTrue(result.app_t0.eq(self.native.app_t0).all())
        self.assertTrue(result.aki_dynamic_first.eq(self.native.aki_dynamic_first).all())
        self.assertTrue(result.hospital_expire_flag.eq(0).all())

    def test_cohort_SQL_requires_source_read_only(self):
        q1=(STAGES/"01_cohort_819_candidates.sql").read_text("utf-8").upper()
        q2=(STAGES/"02_cohort_641_first_adult_stay.sql").read_text("utf-8").upper()
        self.assertIn("223763",q1)
        self.assertIn("INTERVAL '12 HOURS'",q1)
        self.assertIn("DENSE_RANK()",q2)
        self.assertEqual(q2.count("%S::BIGINT[]"),1)
        for q in (q1,q2):
            for forbidden in ("INSERT INTO ","UPDATE MIMICIV","DELETE FROM MIMICIV",
                              "DROP TABLE ","TRUNCATE "):
                self.assertNotIn(forbidden,q)


if __name__=="__main__":
    unittest.main()
