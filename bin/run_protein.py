#!/usr/bin/env python3
"""Cis-MR and colocalisation for one protein against the prepared outcome."""
import argparse
import json

from pdmr.analysis import analyse_protein

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--protein-id", required=True)
ap.add_argument("--gene", required=True)
ap.add_argument("--chr", required=True)
ap.add_argument("--start", type=int, required=True, help="gene start (GRCh38)")
ap.add_argument("--end", type=int, required=True, help="gene end (GRCh38)")
ap.add_argument("--exposure", required=True)
ap.add_argument("--exposure-format", default="ukbppp", choices=["ukbppp", "generic"])
ap.add_argument("--outcome-dir", required=True)
ap.add_argument("--ld-ref", default=None, help="PLINK prefix (GRCh37); may contain {chr}")
ap.add_argument("--prefix", required=True)
for k, t in [("window-kb", float), ("info-min", float), ("maf-min", float), ("p-threshold", float),
             ("f-min", float), ("clump-r2", float), ("clump-kb", float), ("p1", float), ("p2", float),
             ("p12", float), ("outcome-ncase", float), ("outcome-ncontrol", float)]:
    ap.add_argument(f"--{k}", type=t, default=None)
a = ap.parse_args()

meta = {"protein_id": a.protein_id, "gene": a.gene, "chr": a.chr, "start": a.start, "end": a.end}
opts = {k: getattr(a, k) for k in ["window_kb", "info_min", "maf_min", "p_threshold", "f_min", "clump_r2",
                                    "clump_kb", "p1", "p2", "p12", "outcome_ncase", "outcome_ncontrol"]}
ld = a.ld_ref if a.ld_ref and a.ld_ref not in ("NO_LD_REF", "null") else None
summary, mr, region, inst = analyse_protein(meta, a.exposure, a.outcome_dir, ld_ref=ld,
                                            exposure_format=a.exposure_format, **opts)
with open(f"{a.prefix}.summary.json", "w") as fh:
    json.dump(summary, fh, default=float, indent=1)
mr.to_csv(f"{a.prefix}.mr.tsv", sep="\t", index=False)
region.to_csv(f"{a.prefix}.region.tsv.gz", sep="\t", index=False)
inst.to_csv(f"{a.prefix}.instruments.tsv", sep="\t", index=False)
print(f"{a.gene}: {summary['status']}  MR p={summary.get('mr_p')}  PP.H4={summary.get('coloc_PP.H4')}")
