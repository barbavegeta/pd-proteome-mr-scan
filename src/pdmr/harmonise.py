"""Align outcome effects to the exposure effect allele.

Rules follow TwoSampleMR ``harmonise_data(action = 2)``: matching alleles are
kept, swapped alleles have the outcome effect flipped, strand-flipped
non-palindromic variants are complemented, and palindromic (A/T, C/G)
variants are kept only when their frequency makes the strand inferable
(MAF below ``palindrome_maf``, default 0.42) and both studies agree on which
allele is minor.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGT", "TGCA")


def complement(a: pd.Series) -> pd.Series:
    return a.str.translate(COMP)


def is_palindromic(a1: pd.Series, a2: pd.Series) -> pd.Series:
    return (a1.str.len() == 1) & (a2.str.len() == 1) & (complement(a1) == a2)


def harmonise(exp: pd.DataFrame, out: pd.DataFrame, palindrome_maf=0.42):
    """Return (merged table, audit dict). Outcome columns get the ``_out`` suffix."""
    audit = {"n_exposure": int(len(exp)), "n_outcome_window": int(len(out))}
    m = exp.merge(out, on=["chr", "pos37"], suffixes=("_exp", "_out"))
    audit["n_position_match"] = int(len(m))
    if m.empty:
        audit.update(n_harmonised=0, n_palindromic_dropped=0, n_allele_mismatch=0)
        return m, audit

    same = (m.ea_exp == m.ea_out) & (m.oa_exp == m.oa_out)
    swap = (m.ea_exp == m.oa_out) & (m.oa_exp == m.ea_out)
    snv = (m.ea_exp.str.len() == 1) & (m.oa_exp.str.len() == 1)
    c_ea, c_oa = complement(m.ea_out), complement(m.oa_out)
    strand_same = snv & ~same & ~swap & (m.ea_exp == c_ea) & (m.oa_exp == c_oa)
    strand_swap = snv & ~same & ~swap & (m.ea_exp == c_oa) & (m.oa_exp == c_ea)
    ok = same | swap | strand_same | strand_swap
    audit["n_allele_mismatch"] = int((~ok).sum())
    m = m[ok].copy()
    flip = (swap | strand_swap)[ok]
    m.loc[flip, "beta_out"] = -m.loc[flip, "beta_out"]
    m.loc[flip, "eaf_out"] = 1 - m.loc[flip, "eaf_out"]

    pal = is_palindromic(m.ea_exp, m.oa_exp)
    maf_exp = np.minimum(m.eaf_exp, 1 - m.eaf_exp)
    # For palindromes the "same"/"swap" call is ambiguous; trust it only if both
    # studies put the effect allele on the same side of 0.5 and MAF is informative.
    agree = (m.eaf_exp - 0.5) * (m.eaf_out - 0.5) > 0
    pal_keep = pal & (maf_exp < palindrome_maf) & agree.fillna(False)
    drop_pal = pal & ~pal_keep
    audit["n_palindromic_dropped"] = int(drop_pal.sum())
    m = m[~drop_pal]

    # Multi-allelic or duplicated rows: keep the strongest exposure association.
    m = (m.sort_values("log10p", ascending=False) if "log10p" in m else m.sort_values("p_exp")).drop_duplicates(["chr", "pos37", "ea_exp", "oa_exp"])
    m = m.sort_values("pos37").reset_index(drop=True)
    audit["n_harmonised"] = int(len(m))
    return m, audit
