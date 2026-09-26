#!/usr/bin/env python3
"""Download UKB-PPP summary statistics one protein at a time and keep only the cis region.

UKB-PPP (Sun et al. 2023, Nature) distributes one tar per protein on Synapse,
each holding one REGENIE file per chromosome. The full European discovery set
is well over a terabyte, so this script streams it: download a protein's tar,
extract only the chromosome that carries the gene, keep the cis window, write a
small gzipped file, and delete the tar before moving on. Re-running skips
proteins that are already done.

Requirements: ``pip install synapseclient`` and a Synapse account that has
accepted the UKB-PPP terms. Set SYNAPSE_AUTH_TOKEN (personal access token).
Find the folder ID of the European discovery summary statistics on the UKB-PPP
Synapse page and pass it with --folder-id.

    python scripts/fetch_ukbppp_cis.py --folder-id synXXXXXXXX \\
        --protein-map assets/olink_protein_map_3k_v1.tsv --proteins proteins.txt \\
        --outdir data/ukbppp_cis --window-kb 500

NOTE: the Synapse download step could not be exercised in the development
sandbox (no network access to Synapse); the tar extraction and slicing are
covered by tests/test_fetch_extract.py.
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import sys
import tarfile
from pathlib import Path

import pandas as pd


def cis_slice_from_tar(tar_path, chrom, start38, end38, out_path):
    """Extract the chromosome file for ``chrom`` from a UKB-PPP tar and keep GENPOS in [start, end]."""
    with tarfile.open(tar_path) as tf:
        members = [m for m in tf.getmembers()
                   if m.isfile() and f"chr{chrom}_" in Path(m.name).name and m.name.endswith(".gz")]
        if not members:
            raise FileNotFoundError(f"no chr{chrom} file inside {tar_path}")
        raw = tf.extractfile(members[0]).read()
    keep = []
    reader = pd.read_csv(io.BytesIO(raw), sep=r"\s+", compression="gzip", chunksize=500_000,
                         dtype={"CHROM": str, "ID": str, "ALLELE0": str, "ALLELE1": str})
    for chunk in reader:
        keep.append(chunk[chunk["GENPOS"].between(start38, end38)])
    out = pd.concat(keep, ignore_index=True)
    with gzip.open(out_path, "wt") as fh:
        out.to_csv(fh, sep=" ", index=False)
    return len(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--folder-id", required=True, help="Synapse ID of the UKB-PPP European discovery folder")
    ap.add_argument("--protein-map", required=True)
    ap.add_argument("--proteins", help="text file of UKBPPP_ProteinIDs or gene symbols (default: all)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--window-kb", type=float, default=500)
    ap.add_argument("--tmpdir", default=None, help="where tars are downloaded (needs a few GB free)")
    a = ap.parse_args()

    import synapseclient  # imported here so the rest of the module works without it

    syn = synapseclient.Synapse()
    syn.login(authToken=os.environ.get("SYNAPSE_AUTH_TOKEN"))

    pm = pd.read_csv(a.protein_map, sep="\t", dtype={"chr": str})
    if a.proteins:
        wanted = {l.strip() for l in open(a.proteins) if l.strip()}
        pm = pm[pm.UKBPPP_ProteinID.isin(wanted) | pm["HGNC.symbol"].isin(wanted)]
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tmp = Path(a.tmpdir or outdir / "_tmp")
    tmp.mkdir(exist_ok=True)

    children = {c["name"]: c["id"] for c in syn.getChildren(a.folder_id, includeTypes=["file"])}
    w = int(a.window_kb * 1000)
    for r in pm.itertuples():
        pid = r.UKBPPP_ProteinID
        out = outdir / f"{pid.replace(':', '_')}.cis.regenie.gz"
        if out.exists():
            continue
        stem = pid.replace(":", "_")
        match = [n for n in children if n.startswith(stem) and n.endswith(".tar")]
        if not match:
            print(f"[skip] {pid}: no tar named {stem}*.tar in folder", file=sys.stderr)
            continue
        ent = syn.get(children[match[0]], downloadLocation=str(tmp))
        try:
            n = cis_slice_from_tar(ent.path, str(r.chr).replace("X", "23"),
                                   max(0, int(r.gene_start) - w), int(r.gene_end) + w, out)
            print(f"[ok] {pid}: {n} cis variants")
        finally:
            Path(ent.path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
