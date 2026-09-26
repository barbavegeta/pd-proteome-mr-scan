"""Readers that turn pQTL and GWAS summary statistics into one standard table.

Standard columns: chr, pos37, ea, oa, eaf, beta, se, p, n (+ source extras).
``pos37`` (GRCh37) is the join key because the PD GWAS (ieu-b-7) and the
1000 Genomes LD reference are GRCh37, while UKB-PPP reports GRCh38 in GENPOS
and GRCh37 inside the variant ID.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

STD = ["chr", "pos37", "ea", "oa", "eaf", "beta", "se", "p", "n"]


def _open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def norm_chr(x):
    return str(x).replace("chr", "").upper().replace("23", "X")


# --------------------------------------------------------------------------- exposure
def read_ukbppp_regenie(path, chrom, start38, end38, info_min=0.8, maf_min=0.01, chunksize=500_000):
    """Read a UKB-PPP REGENIE file and keep one cis window (GRCh38 coordinates).

    UKB-PPP rows look like ``7 22800054 7:22839673:A:G:imp:v1 A G 0.578 ...``:
    GENPOS is GRCh38, the ID carries the GRCh37 position, BETA and A1FREQ refer
    to ALLELE1.
    """
    usecols = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "LOG10P"]
    keep = []
    chrom = norm_chr(chrom)
    reader = pd.read_csv(path, sep=r"\s+", usecols=usecols, chunksize=chunksize,
                         dtype={"CHROM": str, "ID": str, "ALLELE0": str, "ALLELE1": str})
    for chunk in reader:
        chunk = chunk[(chunk["CHROM"].map(norm_chr) == chrom)
                      & chunk["GENPOS"].between(start38, end38)]
        if len(chunk):
            keep.append(chunk)
    if not keep:
        return pd.DataFrame(columns=STD + ["pos38", "info", "variant_id"])
    d = pd.concat(keep, ignore_index=True)
    parts = d["ID"].str.split(":", expand=True)
    out = pd.DataFrame({
        "chr": chrom,
        "pos37": pd.to_numeric(parts[1], errors="coerce"),
        "pos38": d["GENPOS"].astype(int),
        "ea": d["ALLELE1"].str.upper(),
        "oa": d["ALLELE0"].str.upper(),
        "eaf": d["A1FREQ"].astype(float),
        "beta": d["BETA"].astype(float),
        "se": d["SE"].astype(float),
        "p": np.power(10.0, -d["LOG10P"].astype(float)),
        "log10p": d["LOG10P"].astype(float),
        "n": d["N"].astype(float),
        "info": d["INFO"].astype(float),
        "variant_id": d["ID"],
    })
    out = out.dropna(subset=["pos37", "beta", "se"])
    out["pos37"] = out["pos37"].astype(int)
    maf = np.minimum(out["eaf"], 1 - out["eaf"])
    out = out[(out["info"] >= info_min) & (maf >= maf_min) & (out["se"] > 0)]
    return out.reset_index(drop=True)


def read_generic(path, chrom, start, end, pos_col="pos37", **_):
    """Generic whitespace/tab table with the standard columns (for other pQTL sources)."""
    d = pd.read_csv(path, sep=r"\s+", dtype={"chr": str})
    d["chr"] = d["chr"].map(norm_chr)
    d = d[(d["chr"] == norm_chr(chrom)) & d[pos_col].between(start, end)]
    return d.reset_index(drop=True)


# --------------------------------------------------------------------------- outcome
def iter_opengwas_vcf(path, chunksize=1_000_000):
    """Yield standard-format chunks from an IEU OpenGWAS VCF (ES/SE/LP/AF/SS/NC/ID)."""
    with _open_text(path) as fh:
        n_meta = 0
        for line in fh:
            if line.startswith("##"):
                n_meta += 1
                continue
            header = line.lstrip("#").rstrip("\n").split("\t")
            break
    reader = pd.read_csv(path, sep="\t", skiprows=n_meta + 1, names=header, chunksize=chunksize,
                         dtype={"CHROM": str, "ID": str, "REF": str, "ALT": str}, compression="infer")
    for chunk in reader:
        fmt = chunk["FORMAT"].iloc[0].split(":")
        vals = chunk.iloc[:, -1].str.split(":", expand=True)
        vals.columns = fmt[: vals.shape[1]]
        out = pd.DataFrame({
            "chr": chunk["CHROM"].map(norm_chr),
            "pos37": chunk["POS"].astype(int),
            "rsid": chunk["ID"],
            "ea": chunk["ALT"].str.upper(),
            "oa": chunk["REF"].str.upper(),
            "eaf": pd.to_numeric(vals.get("AF"), errors="coerce"),
            "beta": pd.to_numeric(vals["ES"], errors="coerce"),
            "se": pd.to_numeric(vals["SE"], errors="coerce"),
            "p": np.power(10.0, -pd.to_numeric(vals["LP"], errors="coerce")),
            "n": pd.to_numeric(vals.get("SS"), errors="coerce"),
            "ncase": pd.to_numeric(vals.get("NC"), errors="coerce"),
        })
        yield out.dropna(subset=["beta", "se"])


def prepare_outcome(path, outdir, fmt="opengwas_vcf"):
    """Split an outcome GWAS into one parquet file per chromosome for fast cis lookups."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    buckets: dict[str, list[pd.DataFrame]] = {}
    if fmt == "opengwas_vcf":
        chunks = iter_opengwas_vcf(path)
    elif fmt == "generic":
        chunks = [pd.read_csv(path, sep=r"\s+", dtype={"chr": str})]
    else:
        raise ValueError(f"unknown outcome format {fmt}")
    for chunk in chunks:
        chunk["chr"] = chunk["chr"].map(norm_chr)
        for c, g in chunk.groupby("chr"):
            buckets.setdefault(c, []).append(g)
    written = []
    for c, parts in buckets.items():
        f = outdir / f"chr{c}.parquet"
        pd.concat(parts, ignore_index=True).sort_values("pos37").to_parquet(f, index=False)
        written.append(f)
    return written


def read_outcome_window(outcome_dir, chrom, start37, end37):
    f = Path(outcome_dir) / f"chr{norm_chr(chrom)}.parquet"
    if not f.exists():
        return pd.DataFrame(columns=STD + ["rsid"])
    d = pd.read_parquet(f, filters=[("pos37", ">=", int(start37)), ("pos37", "<=", int(end37))])
    return d.reset_index(drop=True)


READERS = {"ukbppp": read_ukbppp_regenie, "generic": read_generic}
