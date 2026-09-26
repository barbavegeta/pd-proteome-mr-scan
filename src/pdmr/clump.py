"""LD clumping against a PLINK reference panel (e.g. 1000 Genomes EUR, GRCh37).

Greedy clumping as in PLINK ``--clump``: take the most significant remaining
variant, drop everything within ``kb`` that has r^2 above ``r2`` with it, and
repeat. Variants missing from the reference are dropped, matching
``TwoSampleMR::clump_data`` / ``ieugwasr::ld_clump`` behaviour.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


class LDReference:
    def __init__(self, prefix: str):
        self.prefix = str(prefix)
        self._bim = {}

    def _paths(self, chrom):
        p = self.prefix.replace("{chr}", str(chrom))
        return p + ".bed", p + ".bim", p + ".fam"

    def available(self, chrom) -> bool:
        return all(Path(f).exists() for f in self._paths(chrom))

    def bim(self, chrom) -> pd.DataFrame:
        if chrom not in self._bim:
            _, bim, _ = self._paths(chrom)
            b = pd.read_csv(bim, sep=r"\s+", header=None, names=["chr", "snp", "cm", "pos37", "a1", "a2"],
                            dtype={"chr": str, "a1": str, "a2": str})
            b["chr"] = b["chr"].str.replace("chr", "")
            b["idx"] = np.arange(len(b))
            self._bim[chrom] = b[b["chr"] == str(chrom)]
        return self._bim[chrom]

    def genotypes(self, chrom, idx) -> np.ndarray:
        from bed_reader import open_bed
        bed, _, _ = self._paths(chrom)
        with open_bed(bed) as b:
            g = b.read(index=np.s_[:, np.asarray(idx)], dtype="float32")
        col_mean = np.nanmean(g, axis=0)
        g = np.where(np.isnan(g), col_mean, g)
        return g

    def match(self, chrom, variants: pd.DataFrame) -> pd.Series:
        """Return reference column index for each variant (by position and allele pair)."""
        b = self.bim(chrom)
        m = variants[["pos37", "ea", "oa"]].reset_index().merge(b[["pos37", "a1", "a2", "idx"]], on="pos37")
        ok = ((m.ea == m.a1) & (m.oa == m.a2)) | ((m.ea == m.a2) & (m.oa == m.a1))
        m = m[ok].drop_duplicates("index")
        return m.set_index("index")["idx"].reindex(variants.index)


def r2_matrix(g: np.ndarray) -> np.ndarray:
    g = g - g.mean(axis=0)
    sd = g.std(axis=0)
    sd[sd == 0] = np.inf
    z = g / sd
    r = (z.T @ z) / g.shape[0]
    return r ** 2


def clump(cands: pd.DataFrame, ref: LDReference, chrom, r2=0.001, kb=1000):
    """Greedy clump; ``cands`` needs pos37, ea, oa, p. Returns (kept frame, n_missing_from_ref)."""
    idx = ref.match(chrom, cands)
    present = idx.notna()
    c = cands[present].copy()
    c["ref_idx"] = idx[present].astype(int)
    if c.empty:
        return c, int((~present).sum())
    # -log10(p) avoids ties when p underflows to 0 for very strong pQTLs.
    c = (c.sort_values("log10p", ascending=False) if "log10p" in c else c.sort_values("p")).reset_index(drop=True)
    R2 = r2_matrix(ref.genotypes(chrom, c["ref_idx"].values))
    pos = c["pos37"].values
    alive = np.ones(len(c), bool)
    keep = []
    for i in range(len(c)):
        if not alive[i]:
            continue
        keep.append(i)
        near = np.abs(pos - pos[i]) <= kb * 1000
        alive &= ~(near & (R2[i] > r2))
        alive[i] = False
    return c.iloc[keep].reset_index(drop=True), int((~present).sum())
