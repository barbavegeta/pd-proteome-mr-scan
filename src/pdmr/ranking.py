"""Combine per-protein results, control FDR, and assign evidence tiers."""
from __future__ import annotations

import numpy as np
import pandas as pd

TIERS = [
    "1: MR + colocalisation",
    "2: MR + suggestive colocalisation",
    "3: MR, distinct causal variants",
    "4: MR, colocalisation inconclusive",
    "No MR evidence",
    "Not testable",
]


def bh_fdr(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q-values; NaN p-values are ignored."""
    q = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if not ok.any():
        return q
    pv = p[ok].values
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def assign_tier(row, fdr=0.05, h4_strong=0.8, h4_suggestive=0.5, h3_distinct=0.5):
    if pd.isna(row.get("mr_p")):
        return TIERS[5]
    if row["mr_q"] >= fdr:
        return TIERS[4]
    h4, h3 = row.get("coloc_PP.H4", np.nan), row.get("coloc_PP.H3", np.nan)
    if h4 >= h4_strong:
        return TIERS[0]
    if h4 >= h4_suggestive:
        return TIERS[1]
    if h3 >= h3_distinct:
        return TIERS[2]
    return TIERS[3]


def rank(summary: pd.DataFrame, fdr=0.05, **kw) -> pd.DataFrame:
    d = summary.copy()
    d["mr_q"] = bh_fdr(d["mr_p"]) if "mr_p" in d else np.nan
    d["tier"] = d.apply(assign_tier, axis=1, fdr=fdr, **kw)
    d["tier_rank"] = d["tier"].map({t: i for i, t in enumerate(TIERS)})
    d = d.sort_values(["tier_rank", "coloc_PP.H4", "mr_p"], ascending=[True, False, True])
    return d.drop(columns="tier_rank").reset_index(drop=True)
