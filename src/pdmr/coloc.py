"""Approximate Bayes factor colocalisation (Giambartolomei et al. 2014).

A NumPy re-implementation of ``coloc::coloc.abf`` for the single-causal-variant
model. It reproduces the R package's defaults:

* quantitative traits use a prior effect-size SD of ``0.15 * sdY``; when
  ``sdY`` is not supplied it is estimated from ``varbeta``, ``MAF`` and ``N``
  exactly as ``coloc:::sdY.est`` does (regression of 2*N*MAF*(1-MAF) on
  1/varbeta through the origin);
* case-control traits use a prior SD of 0.2 on the log-odds scale;
* posterior probabilities are combined as in ``coloc:::combine.abf``.

Validated against the R package in ``tests/test_coloc.py`` using the
dissertation's SomaScan and UKB-PPP Olink GPNMB inputs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

HYPOTHESES = ("H0", "H1", "H2", "H3", "H4")


@dataclass
class Dataset:
    beta: np.ndarray
    varbeta: np.ndarray
    type: str  # "quant" or "cc"
    N: float | None = None
    MAF: np.ndarray | None = None
    sdY: float | None = None
    s: float | None = None  # case fraction, informational for "cc"

    def prior_sd(self) -> float:
        if self.type == "cc":
            return 0.2
        if self.type != "quant":
            raise ValueError(f"unknown dataset type {self.type!r}")
        sdy = self.sdY if self.sdY is not None else estimate_sdY(self.varbeta, self.MAF, self.N)
        return 0.15 * sdy


def estimate_sdY(varbeta: np.ndarray, maf: np.ndarray | None, n: float | None) -> float:
    """Port of ``coloc:::sdY.est``: sqrt of the slope of nvx ~ 1/varbeta - 1."""
    if maf is None or n is None:
        raise ValueError("sdY not given: MAF and N are required to estimate it")
    varbeta = np.asarray(varbeta, float)
    maf = np.asarray(maf, float)
    oneover = 1.0 / varbeta
    nvx = 2.0 * n * maf * (1.0 - maf)
    slope = np.sum(oneover * nvx) / np.sum(oneover ** 2)
    if slope < 0:
        raise ValueError("estimated sdY is negative; supply sdY explicitly")
    return float(np.sqrt(slope))


def log_abf(beta: np.ndarray, varbeta: np.ndarray, prior_sd: float) -> np.ndarray:
    """Wakefield log approximate Bayes factor for each variant."""
    beta = np.asarray(beta, float)
    v = np.asarray(varbeta, float)
    w2 = prior_sd ** 2
    r = w2 / (w2 + v)
    z = beta / np.sqrt(v)
    return 0.5 * (np.log(1.0 - r) + r * z ** 2)


def _logdiff(a: float, b: float) -> float:
    """log(exp(a) - exp(b)) for a > b, as in coloc:::logdiff."""
    mx = max(a, b)
    return mx + np.log(np.exp(a - mx) - np.exp(b - mx))


def combine_abf(l1: np.ndarray, l2: np.ndarray, p1=1e-4, p2=1e-4, p12=1e-5) -> dict:
    lsum = l1 + l2
    lh0 = 0.0
    lh1 = np.log(p1) + logsumexp(l1)
    lh2 = np.log(p2) + logsumexp(l2)
    lh3 = np.log(p1) + np.log(p2) + _logdiff(logsumexp(l1) + logsumexp(l2), logsumexp(lsum))
    lh4 = np.log(p12) + logsumexp(lsum)
    all_abf = np.array([lh0, lh1, lh2, lh3, lh4])
    pp = np.exp(all_abf - logsumexp(all_abf))
    return {f"PP.{h}": float(v) for h, v in zip(HYPOTHESES, pp)}


def coloc_abf(d1: Dataset, d2: Dataset, p1=1e-4, p2=1e-4, p12=1e-5, return_snp=False):
    """Run colocalisation on two aligned datasets (same variants, same order).

    Returns a dict with ``nsnps`` and ``PP.H0``..``PP.H4``; with
    ``return_snp=True`` also returns per-variant ``SNP.PP.H4``.
    """
    l1 = log_abf(d1.beta, d1.varbeta, d1.prior_sd())
    l2 = log_abf(d2.beta, d2.varbeta, d2.prior_sd())
    if len(l1) != len(l2):
        raise ValueError("datasets must be aligned to the same variants")
    res = {"nsnps": int(len(l1)), **combine_abf(l1, l2, p1, p2, p12)}
    if return_snp:
        lsum = l1 + l2
        return res, np.exp(lsum - logsumexp(lsum))
    return res


def prior_sensitivity(d1: Dataset, d2: Dataset, p12_grid=(1e-6, 5e-6, 1e-5, 5e-5), p1=1e-4, p2=1e-4):
    """PP.H4 across a grid of p12 priors (the posterior for H3/H4 is prior-sensitive)."""
    l1 = log_abf(d1.beta, d1.varbeta, d1.prior_sd())
    l2 = log_abf(d2.beta, d2.varbeta, d2.prior_sd())
    return {float(p): combine_abf(l1, l2, p1, p2, p)["PP.H4"] for p in p12_grid}
