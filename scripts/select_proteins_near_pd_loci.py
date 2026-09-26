#!/usr/bin/env python3
"""Pick the proteins whose genes lie near genome-wide significant PD loci.

A practical first pass before the full proteome: it needs far less data to
download. Loci are defined from the prepared outcome (per-chromosome parquet
from prepare_outcome.py) by merging significant variants within --merge-kb.
The protein map is GRCh38 and ieu-b-7 is GRCh37, so the default margin is
deliberately generous (1.5 Mb); the scan itself matches variants exactly.

    python scripts/select_proteins_near_pd_loci.py --outcome-dir outcome_by_chr \\
        --protein-map assets/olink_protein_map_3k_v1.tsv --out proteins_pd_loci.txt
"""
import argparse
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--outcome-dir", required=True)
ap.add_argument("--protein-map", required=True)
ap.add_argument("--p", type=float, default=5e-8)
ap.add_argument("--merge-kb", type=float, default=500)
ap.add_argument("--margin-kb", type=float, default=1500)
ap.add_argument("--out", default="proteins_pd_loci.txt")
a = ap.parse_args()

loci = []
for f in sorted(Path(a.outcome_dir).glob("chr*.parquet")):
    d = pd.read_parquet(f, columns=["chr", "pos37", "p"])
    d = d[d.p < a.p].sort_values("pos37")
    if d.empty:
        continue
    start = prev = d.pos37.iloc[0]
    for pos in d.pos37.iloc[1:]:
        if pos - prev > a.merge_kb * 1000:
            loci.append((d.chr.iloc[0], start, prev))
            start = pos
        prev = pos
    loci.append((d.chr.iloc[0], start, prev))
print(f"{len(loci)} PD loci at p < {a.p:g}")

pm = pd.read_csv(a.protein_map, sep="\t", dtype={"chr": str})
m = a.margin_kb * 1000
keep = set()
for c, s, e in loci:
    near = pm[(pm.chr == str(c)) & (pm.gene_end >= s - m) & (pm.gene_start <= e + m)]
    keep.update(near.UKBPPP_ProteinID)
Path(a.out).write_text("\n".join(sorted(keep)) + "\n")
print(f"{len(keep)} proteins written to {a.out}")
