#!/usr/bin/env python3
"""Render the HTML report (figures embedded) from the aggregated scan results."""
import argparse

from pdmr.report import build_report

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--results", required=True)
ap.add_argument("--region-dir", default=".", help="folder with <gene>.region.tsv.gz / .instruments.tsv")
ap.add_argument("--highlight", default="GPNMB", help="comma-separated genes always shown")
ap.add_argument("--params", default=None, help="JSON with run parameters to print in the report")
ap.add_argument("--out", default="pd_proteome_mr_report.html")
a = ap.parse_args()
build_report(a.results, a.region_dir, a.out, [g for g in a.highlight.split(",") if g], a.params)
print(a.out)
