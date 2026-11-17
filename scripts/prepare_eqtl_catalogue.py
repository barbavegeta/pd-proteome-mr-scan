#!/usr/bin/env python3
"""Split an eQTL Catalogue dataset into one cis slice per molecular trait, for --exposure_format eqtlcat.

Built for the Sun et al. 2018 SomaScan plasma pQTLs (eQTL Catalogue QTD000584,
3,301 INTERVAL participants), which, unlike UKB-PPP, can be downloaded without a
login. The eQTL Catalogue file is GRCh38 and holds every molecular trait's cis
window, sorted by position. This script streams it once and:

  * looks up each variant's GRCh37 position through its rsID in the prepared
    outcome GWAS (the pipeline matches exposure and outcome on GRCh37 position),
    dropping variants the outcome does not have;
  * writes ALT-allele frequency (ac / an), N (an / 2), imputation r2 and an exact
    -log10 p computed from beta / se, which stays finite for very strong pQTLs;
  * writes one TSV per molecular trait;
  * fetches gene coordinates (GRCh38) from the Ensembl REST API and writes the
    samplesheet.

    python scripts/prepare_eqtl_catalogue.py \\
        --sumstats QTD000584.all.tsv.gz --outcome-dir outcome_by_chr \\
        --outdir data/sun2018_cis --samplesheet samplesheet_sun2018.csv
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

COLS = ["molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "pvalue",
        "beta", "se", "ac", "an", "r2", "gene_id", "rsid"]
SLICE_COLS = ["chr", "pos38", "pos37", "rsid", "variant_id", "ea", "oa", "eaf", "beta", "se",
              "p", "log10p", "n", "info"]


def seq_id(trait_id: str) -> str:
    """SomaScan SeqId from an eQTL Catalogue aptamer ID.

    'SNCA.8458.111.3..1' -> '8458-111'; 'AKT1.AKT2.AKT3.3392.68.2..3' -> '3392-68'.
    """
    m = re.search(r"\.(\d+)\.(\d+)\.\d+\.\.\d+$", trait_id)
    return f"{m.group(1)}-{m.group(2)}" if m else trait_id


def safe_name(trait_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", trait_id)


def rsid_to_pos37(outcome_dir: Path, chrom: str) -> pd.Series:
    f = outcome_dir / f"chr{chrom}.parquet"
    if not f.exists():
        return pd.Series(dtype="int64")
    o = pd.read_parquet(f, columns=["rsid", "pos37"])
    o = o[o.rsid.str.startswith("rs", na=False)].drop_duplicates("rsid")
    return o.set_index("rsid")["pos37"]


def ensembl_lookup(gene_ids: list[str], batch=1000) -> pd.DataFrame:
    rows = []
    for i in range(0, len(gene_ids), batch):
        body = json.dumps({"ids": gene_ids[i:i + batch]}).encode()
        req = urllib.request.Request("https://rest.ensembl.org/lookup/id", data=body,
                                     headers={"Content-Type": "application/json",
                                              "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            res = json.load(r)
        for gid, g in res.items():
            if g:
                rows.append({"gene_id": gid, "symbol": g.get("display_name") or gid,
                             "chr": str(g["seq_region_name"]), "start": int(g["start"]),
                             "end": int(g["end"])})
        time.sleep(1)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sumstats", required=True, help="eQTL Catalogue *.all.tsv.gz")
    ap.add_argument("--outcome-dir", required=True, help="output of bin/prepare_outcome.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--chunksize", type=int, default=2_000_000)
    a = ap.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outcome_dir = Path(a.outcome_dir)
    maps: dict[str, pd.Series] = {}
    started: set[str] = set()
    trait_gene: dict[str, str] = {}
    n_in = n_kept = 0

    reader = pd.read_csv(a.sumstats, sep="\t", usecols=COLS, chunksize=a.chunksize,
                         dtype={"chromosome": str, "rsid": str, "ref": str, "alt": str})
    for chunk in reader:
        n_in += len(chunk)
        trait_gene.update(chunk.drop_duplicates("molecular_trait_id")
                          .set_index("molecular_trait_id")["gene_id"].to_dict())
        chunk = chunk.dropna(subset=["rsid", "beta", "se", "an"])
        chunk = chunk[(chunk.se > 0) & (chunk.an > 0)]
        parts = []
        for chrom, g in chunk.groupby("chromosome"):
            chrom = chrom.replace("chr", "")
            if chrom not in maps:
                maps[chrom] = rsid_to_pos37(outcome_dir, chrom)
            pos37 = g.rsid.map(maps[chrom])
            g = g.assign(pos37=pos37, chromosome=chrom).dropna(subset=["pos37"])
            parts.append(g)
        if not parts:
            continue
        d = pd.concat(parts)
        z = (d.beta / d.se).abs()
        log10p = -(norm.logsf(z) + np.log(2)) / np.log(10)
        out = pd.DataFrame({
            "chr": d.chromosome, "pos38": d.position.astype(int), "pos37": d.pos37.astype(int),
            "rsid": d.rsid, "variant_id": d.variant, "ea": d.alt.str.upper(), "oa": d.ref.str.upper(),
            "eaf": d.ac / d.an, "beta": d.beta, "se": d.se, "p": np.power(10.0, -log10p),
            "log10p": log10p, "n": d.an / 2, "info": d.r2.fillna(1.0),
        })
        out["trait"] = d.molecular_trait_id.values
        n_kept += len(out)
        for trait, g in out.groupby("trait"):
            f = outdir / f"{safe_name(trait)}.tsv"
            g[SLICE_COLS].to_csv(f, sep="\t", index=False, mode="a", header=trait not in started)
            started.add(trait)
        print(f"  {n_in:,} rows read, {n_kept:,} kept, {len(started):,} traits", flush=True)

    # Gene coordinates (GRCh38) from Ensembl.
    tg = pd.DataFrame({"protein_id": list(started)})
    tg["gene_id"] = tg.protein_id.map(trait_gene)
    genes = ensembl_lookup(sorted(tg.gene_id.dropna().unique()))
    sheet = tg.merge(genes, on="gene_id", how="inner")
    sheet = sheet[sheet.chr.isin([str(c) for c in range(1, 23)])]
    # Several aptamers can target one gene: keep gene labels unique for display.
    dup = sheet.symbol.duplicated(keep=False)
    sheet["gene"] = sheet.symbol
    sheet.loc[dup, "gene"] = [f"{s}_{seq_id(p)}" for s, p in zip(sheet.symbol[dup], sheet.protein_id[dup])]
    # The same aptamer can be listed twice for one gene ('..1', '..3'): add the entry number.
    again = sheet.gene.duplicated(keep=False)
    sheet.loc[again, "gene"] = [f"{g}.{p.rsplit('..', 1)[-1]}" for g, p in zip(sheet.gene[again], sheet.protein_id[again])]
    sheet["exposure"] = [str((outdir / f"{safe_name(t)}.tsv").resolve()) for t in sheet.protein_id]
    sheet = sheet.assign(chr_num=sheet.chr.astype(int)).sort_values(["chr_num", "start"])
    sheet[["protein_id", "gene", "chr", "start", "end", "exposure", "gene_id"]].to_csv(a.samplesheet, index=False)
    print(f"{len(sheet):,} traits with autosomal gene coordinates written to {a.samplesheet}")


if __name__ == "__main__":
    main()
