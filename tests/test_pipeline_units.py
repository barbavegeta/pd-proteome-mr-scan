"""Unit and end-to-end tests for harmonisation, clumping, ranking and the per-protein run."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pdmr.analysis import analyse_protein
from pdmr.clump import LDReference, clump
from pdmr.harmonise import harmonise
from pdmr.io import prepare_outcome
from pdmr.mr import ivw, mr_egger, weighted_median
from pdmr.ranking import TIERS, bh_fdr, rank

ROOT = Path(__file__).parent
DATA = ROOT / "data"
LD = str(DATA / "ld_ref" / "synthetic_eur_chr22")


def _row(pos, ea, oa, eaf, beta, se=0.01, p=1e-3):
    return dict(chr="1", pos37=pos, ea=ea, oa=oa, eaf=eaf, beta=beta, se=se, p=p, n=1000)


def test_harmonise_flips_and_palindromes():
    exp = pd.DataFrame([_row(1, "A", "G", 0.3, 0.1), _row(2, "A", "G", 0.3, 0.1), _row(3, "A", "G", 0.3, 0.1),
                        _row(4, "A", "T", 0.2, 0.1), _row(5, "A", "T", 0.45, 0.1), _row(6, "A", "C", 0.3, 0.1)])
    out = pd.DataFrame([_row(1, "A", "G", 0.31, 0.5), _row(2, "G", "A", 0.69, 0.5), _row(3, "T", "C", 0.3, 0.5),
                        _row(4, "A", "T", 0.21, 0.5), _row(5, "A", "T", 0.46, 0.5), _row(6, "A", "G", 0.3, 0.5)])
    h, audit = harmonise(exp, out)
    got = dict(zip(h.pos37, h.beta_out))
    assert got[1] == 0.5          # same alleles
    assert got[2] == -0.5         # swapped -> flipped
    assert got[3] == 0.5          # strand flip (T/C == complement of A/G)
    assert got[4] == 0.5          # palindromic but low MAF and frequencies agree -> kept
    assert 5 not in got           # palindromic with MAF 0.45 -> dropped
    assert 6 not in got           # allele mismatch
    assert audit["n_palindromic_dropped"] == 1 and audit["n_allele_mismatch"] == 1


def test_bh_fdr_matches_reference():
    p = pd.Series([0.01, 0.04, 0.03, 0.2, np.nan])
    q = bh_fdr(p)
    assert q.iloc[:4].round(6).tolist() == [0.04, 0.053333, 0.053333, 0.2]
    assert np.isnan(q.iloc[4])


def test_mr_recovers_simulated_causal_effect():
    rng = np.random.default_rng(0)
    bx = rng.uniform(0.05, 0.2, 20)
    sx = np.full(20, 0.005)
    by = 0.3 * bx + rng.normal(0, 0.005, 20)
    sy = np.full(20, 0.005)
    fe, re = ivw(bx, sx, by, sy)
    assert fe["beta"] == pytest.approx(0.3, abs=0.03)
    assert mr_egger(bx, sx, by, sy)["intercept_p"] > 0.01
    assert weighted_median(bx, sx, by, sy)["beta"] == pytest.approx(0.3, abs=0.05)


def test_clumping_keeps_independent_signals():
    ref = LDReference(LD)
    bim = ref.bim("22")
    # take a dense run of variants from one simulated locus: clumping must thin them
    c = bim.iloc[:60].rename(columns={"a1": "ea", "a2": "oa"})[["pos37", "ea", "oa"]].copy()
    c["p"] = np.linspace(1e-20, 1e-9, len(c))
    kept, missing = clump(c, ref, "22", r2=0.01, kb=1000)
    assert missing == 0
    assert 1 <= len(kept) < len(c)
    g = ref.genotypes("22", kept.ref_idx.values)
    r2 = np.corrcoef(g.T) ** 2
    np.fill_diagonal(r2, 0)
    assert r2.max() <= 0.01 + 1e-6


@pytest.fixture(scope="module")
def outcome_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("outcome")
    prepare_outcome(DATA / "pd_outcome_test.vcf.gz", d)
    return d


@pytest.fixture(scope="module")
def scan(outcome_dir):
    sheet = pd.read_csv(DATA / "samplesheet_test.csv")
    rows = []
    for r in sheet.itertuples():
        meta = r._asdict()
        s, *_ = analyse_protein(meta, DATA / r.exposure, outcome_dir, ld_ref=LD)
        rows.append(s)
    return rank(pd.DataFrame(rows)).set_index("gene")


def test_simulated_scenarios_land_in_expected_tiers(scan):
    assert scan.loc["SYN_SHARED", "tier"] == TIERS[0]
    assert scan.loc["SYN_MULTI", "tier"] == TIERS[0]
    assert scan.loc["SYN_MULTI", "n_instruments"] == 3
    assert scan.loc["SYN_DISTINCT", "coloc_PP.H3"] > 0.9
    assert scan.loc["SYN_DISTINCT", "tier"] in (TIERS[2], TIERS[4])
    assert scan.loc["SYN_NULL", "tier"] == TIERS[4]
    assert scan.loc["SYN_NULL", "coloc_PP.H1"] > 0.8


def test_real_gpnmb_locus_reproduces_dissertation(scan):
    g = scan.loc["GPNMB"]
    # no chr7 LD panel in the test data -> lead-variant Wald ratio, as in the dissertation
    assert g["exposure_lead_variant"] == "7:23306141:G:A:imp:v1"   # rs75801644
    assert g["mr_beta"] == pytest.approx(-0.0671699193245443, rel=1e-9)
    assert g["mr_se"] == pytest.approx(0.04829794465413305, rel=1e-9)
    # distinct signals, as in the dissertation's Olink coloc (PP.H3 = 0.9998)
    assert g["coloc_PP.H3"] == pytest.approx(0.9998, abs=5e-4)
