"""Golden tests: the Python implementation must reproduce the dissertation's R results.

Expected values come from files produced with R 4.5.1, coloc and TwoSampleMR
(see tests/fixtures/dissertation/ and the GPNMB-Parkinson-MR-Colocalisation repository).
"""
from pathlib import Path

import pandas as pd
import pytest

from pdmr.coloc import Dataset, coloc_abf, prior_sensitivity
from pdmr.mr import ivw, wald_ratio

FIX = Path(__file__).parent / "fixtures" / "dissertation"
N_CASE, N_CTRL = 33674, 449056
N_PD = N_CASE + N_CTRL
PP = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]


def test_coloc_somascan_matches_r():
    d = pd.read_csv(FIX / "gpnmb_coloc_merged_harmonised_peer_review.csv")
    d1 = Dataset(d.beta_exposure.values, d.varbeta_exposure.values, "quant", N=996, MAF=d.maf.values)
    d2 = Dataset(d.beta_outcome.values, d.varbeta_outcome.values, "cc", N=N_PD, s=N_CASE / N_PD)
    got = coloc_abf(d1, d2, p12=1e-5)
    exp = pd.read_csv(FIX / "gpnmb_final_coloc_results_peer_review.csv").iloc[0]
    assert got["nsnps"] == exp["nsnps"] == 150
    for k in PP:
        assert got[k] == pytest.approx(exp[f"{k}.abf"], rel=1e-8, abs=1e-15)


def test_coloc_olink_dense_matches_r():
    d = pd.read_csv(FIX / "01_olink_coloc_dense_input_after_technical_qc.csv")
    d1 = Dataset(d.BETA.values, d.varbeta_exposure.values, "quant", N=32743, MAF=d.exposure_MAF.values)
    d2 = Dataset(d.out_beta_to_expA1.values, d.varbeta_outcome.values, "cc", N=N_PD, s=N_CASE / N_PD)
    got = coloc_abf(d1, d2, p12=1e-5)
    exp = pd.read_csv(FIX / "02_olink_coloc_abf_main_summary.csv").iloc[0]
    assert got["nsnps"] == exp["nsnps"] == 4018
    for k in PP:
        assert got[k] == pytest.approx(exp[k.replace(".", "_")], rel=1e-8, abs=1e-15)


def test_coloc_olink_prior_sensitivity_matches_r():
    d = pd.read_csv(FIX / "01_olink_coloc_dense_input_after_technical_qc.csv")
    d1 = Dataset(d.BETA.values, d.varbeta_exposure.values, "quant", N=32743, MAF=d.exposure_MAF.values)
    d2 = Dataset(d.out_beta_to_expA1.values, d.varbeta_outcome.values, "cc", N=N_PD, s=N_CASE / N_PD)
    exp = pd.read_csv(FIX / "04_olink_coloc_abf_prior_sensitivity.csv")
    got = prior_sensitivity(d1, d2, tuple(exp.p12))
    for p12, h4 in zip(exp.p12, exp.PP_H4):
        assert got[float(p12)] == pytest.approx(h4, rel=1e-8)


@pytest.mark.parametrize("setting", ["strict_r2_0p001_10Mb", "standard_r2_0p01_10Mb", "liberal_r2_0p05_1Mb"])
def test_ivw_matches_twosamplemr(setting):
    inst = pd.read_csv(FIX / f"02_clumped_instruments_{setting}.csv")
    fe, re = ivw(inst.exposure_beta, inst.exposure_se, inst.outcome_beta, inst.outcome_se)
    exp = pd.read_csv(FIX / "05_multiplicative_random_effects_ivw_results.csv").set_index("setting").loc[setting]
    assert fe["n_snp"] == exp.n_instruments
    assert fe["beta"] == pytest.approx(exp.fixed_effect_beta, rel=1e-10)
    assert fe["se"] == pytest.approx(exp.fixed_effect_se, rel=1e-10)
    assert fe["p"] == pytest.approx(exp.fixed_effect_p, rel=1e-8)
    assert fe["Q"] == pytest.approx(exp.cochran_Q, rel=1e-10)
    assert fe["Q_p"] == pytest.approx(exp.cochran_Q_p, rel=1e-6)
    assert re["phi"] == pytest.approx(exp.multiplicative_RE_phi, rel=1e-10)
    assert re["se"] == pytest.approx(exp.multiplicative_RE_se, rel=1e-10)
    assert re["p"] == pytest.approx(exp.multiplicative_RE_p, rel=1e-8)


def test_wald_ratio_matches_r():
    r = pd.read_csv(FIX / "01_lead_variant_wald_mr.csv").iloc[0]
    got = wald_ratio(r.exposure_beta, r.exposure_se, r.outcome_beta, r.outcome_se, second_order=True)
    assert got["beta"] == pytest.approx(r.mr_beta, rel=1e-12)
    assert got["se"] == pytest.approx(r.mr_se, rel=1e-12)
    assert got["p"] == pytest.approx(r.mr_p, rel=1e-8)
