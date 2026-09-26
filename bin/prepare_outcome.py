#!/usr/bin/env python3
"""Split the outcome GWAS into per-chromosome parquet files for fast cis-window lookups."""
import argparse

from pdmr.io import prepare_outcome

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--gwas", required=True, help="Outcome GWAS (OpenGWAS VCF, .vcf or .vcf.gz, or generic table)")
ap.add_argument("--format", default="opengwas_vcf", choices=["opengwas_vcf", "generic"])
ap.add_argument("--outdir", default="outcome_by_chr")
a = ap.parse_args()
for f in prepare_outcome(a.gwas, a.outdir, a.format):
    print(f)
