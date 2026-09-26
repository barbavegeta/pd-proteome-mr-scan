#!/usr/bin/env python3
"""Build the pipeline samplesheet from the Olink protein map and a folder of pQTL files.

Each protein is matched to a file in --sumstats-dir whose name contains its
UKBPPP_ProteinID (with ':' or '_' separators) and, when several chromosome
files exist, the chromosome carrying the gene (``chr<N>_``). Works for the
cis slices written by fetch_ukbppp_cis.py and for fully extracted UKB-PPP tars.

    python scripts/build_samplesheet.py --protein-map assets/olink_protein_map_3k_v1.tsv \\
        --sumstats-dir data/ukbppp_cis --out samplesheet.csv
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--protein-map", required=True)
    ap.add_argument("--sumstats-dir", required=True)
    ap.add_argument("--proteins", help="optional text file restricting to these protein IDs or gene symbols")
    ap.add_argument("--out", default="samplesheet.csv")
    a = ap.parse_args()

    pm = pd.read_csv(a.protein_map, sep="\t", dtype={"chr": str})
    if a.proteins:
        wanted = {l.strip() for l in open(a.proteins) if l.strip()}
        pm = pm[pm.UKBPPP_ProteinID.isin(wanted) | pm["HGNC.symbol"].isin(wanted)]
    files = [p for p in Path(a.sumstats_dir).rglob("*") if p.is_file()]
    rows, missing = [], []
    for _, r in pm.iterrows():
        pid = r["UKBPPP_ProteinID"]
        chrom = str(r["chr"]).replace("X", "23")
        hits = [f for f in files if pid in f.name or pid.replace(":", "_") in f.name]
        if len(hits) > 1:
            hits = [f for f in hits if f"chr{chrom}_" in f.name] or hits
        if not hits:
            missing.append(pid)
            continue
        rows.append({"protein_id": pid, "gene": r["HGNC.symbol"], "chr": chrom,
                     "start": int(r["gene_start"]), "end": int(r["gene_end"]),
                     "exposure": str(hits[0].resolve())})
    sheet = pd.DataFrame(rows, columns=["protein_id", "gene", "chr", "start", "end", "exposure"])
    # UKB-PPP measures a few proteins on more than one panel: keep gene names unique for display
    dup = sheet.gene.duplicated(keep=False)
    if dup.any():
        sheet.loc[dup, "gene"] = [f"{g}_{p.split(':')[2]}" for g, p in zip(sheet.gene[dup], sheet.protein_id[dup])]
    sheet.to_csv(a.out, index=False)
    print(f"{len(sheet)} proteins written to {a.out}; {len(missing)} without a matching file")


if __name__ == "__main__":
    main()
