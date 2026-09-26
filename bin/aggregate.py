#!/usr/bin/env python3
"""Merge per-protein results, apply Benjamini-Hochberg FDR and assign evidence tiers."""
import argparse
import glob
import json

import pandas as pd

from pdmr.ranking import rank

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--summaries", nargs="+", required=True)
ap.add_argument("--mr", nargs="*", default=[])
ap.add_argument("--fdr", type=float, default=0.05)
ap.add_argument("--out", default="scan_results.tsv")
a = ap.parse_args()

files = [f for pat in a.summaries for f in glob.glob(pat)]
rows = [json.load(open(f)) for f in files]
res = rank(pd.DataFrame(rows), fdr=a.fdr)
lead = ["tier", "gene", "protein_id", "mr_OR", "mr_OR_lci", "mr_OR_uci", "mr_p", "mr_q", "mr_method",
        "n_instruments", "coloc_PP.H4", "coloc_PP.H3", "coloc_nsnps", "flags", "status"]
res = res[[c for c in lead if c in res] + [c for c in res if c not in lead]]
res.to_csv(a.out, sep="\t", index=False)
mr_files = [f for pat in a.mr for f in glob.glob(pat)]
mr_tabs = [pd.read_csv(f, sep="\t") for f in mr_files if open(f).read().strip()]
if mr_tabs:
    pd.concat(mr_tabs, ignore_index=True).to_csv("mr_all_methods.tsv", sep="\t", index=False)
print(res["tier"].value_counts().to_string())
