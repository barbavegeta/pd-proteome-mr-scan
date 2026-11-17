#!/usr/bin/env python3
"""Download UKB-PPP summary statistics one protein at a time and keep only the cis region.

UKB-PPP (Sun et al. 2023, Nature) distributes one tar per protein on Synapse
(~545 MB each, ~1.6 TB for the European discovery set), each holding one
REGENIE file per chromosome. Only the chromosome carrying the gene is needed.

By default this script never downloads a whole tar. The files sit in S3, which
serves byte ranges, and a tar is a sequence of 512-byte headers followed by
member data. So it reads the headers one by one (a few kB), downloads only the
cis chromosome (~30 MB), keeps the cis window and writes a small gzipped file.
The European discovery set comes down as under 100 GB instead of ~1.6 TB.
``--full-download`` falls back to downloading each tar. Re-running skips
proteins that are already done.

Requirements: ``pip install synapseclient`` and a Synapse account that has
accepted the UKB-PPP terms. Set SYNAPSE_AUTH_TOKEN (personal access token).
The European discovery folder is syn51365303.

    python scripts/fetch_ukbppp_cis.py --folder-id syn51365303 \\
        --protein-map assets/olink_protein_map_3k_v1.tsv --proteins proteins.txt \\
        --outdir data/ukbppp_cis --window-kb 500 --workers 8
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import re
import tarfile
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


def slice_regenie_gz(raw: bytes, start38, end38, out_path):
    """Keep GENPOS in [start, end] from a gzipped REGENIE file and write it gzipped.

    Writes to a temporary file first, so an interrupted run never leaves a
    partial output that a re-run would mistake for a finished protein.
    """
    keep = []
    reader = pd.read_csv(io.BytesIO(raw), sep=r"\s+", compression="gzip", chunksize=500_000,
                         dtype={"CHROM": str, "ID": str, "ALLELE0": str, "ALLELE1": str})
    for chunk in reader:
        keep.append(chunk[chunk["GENPOS"].between(start38, end38)])
    out = pd.concat(keep, ignore_index=True)
    out_path = Path(out_path)
    tmp = out_path.with_name(out_path.name + ".part")
    with gzip.open(tmp, "wt") as fh:
        out.to_csv(fh, sep=" ", index=False)
    tmp.replace(out_path)
    return len(out)


def _is_chrom_file(name, chrom):
    return f"chr{chrom}_" in Path(name).name and name.endswith(".gz")


def cis_slice_from_tar(tar_path, chrom, start38, end38, out_path):
    """Extract the chromosome file for ``chrom`` from a local UKB-PPP tar and keep GENPOS in [start, end]."""
    with tarfile.open(tar_path) as tf:
        members = [m for m in tf.getmembers() if m.isfile() and _is_chrom_file(m.name, chrom)]
        if not members:
            raise FileNotFoundError(f"no chr{chrom} file inside {tar_path}")
        raw = tf.extractfile(members[0]).read()
    return slice_regenie_gz(raw, start38, end38, out_path)


def tar_members_by_range(read_range):
    """List (name, size, data_offset) for a tar's regular files using only byte-range reads.

    ``read_range(start, length)`` returns bytes. One 512-byte header is read per
    member and its data is skipped, so this costs a few kB however large the tar is.
    GNU long names and pax ``path`` records are honoured.
    """
    members, off, pending_name = [], 0, None
    while True:
        hdr = read_range(off, 512)
        if len(hdr) < 512 or hdr == b"\0" * 512:
            return members
        ti = tarfile.TarInfo.frombuf(hdr, tarfile.ENCODING, "surrogateescape")
        data_off = off + 512
        if ti.type == tarfile.GNUTYPE_LONGNAME:
            pending_name = read_range(data_off, ti.size).rstrip(b"\0").decode("utf-8", "replace")
        elif ti.type in (tarfile.XHDTYPE, tarfile.XGLTYPE):
            m = re.search(r"\d+ path=([^\n]*)\n", read_range(data_off, ti.size).decode("utf-8", "replace"))
            if m:
                pending_name = m.group(1)
        else:
            if ti.isfile():
                members.append((pending_name or ti.name, ti.size, data_off))
            pending_name = None
        off = data_off + -(-ti.size // 512) * 512


def cis_slice_by_range(read_range, chrom, start38, end38, out_path):
    """Like cis_slice_from_tar, but reads only the needed member of a remote tar."""
    hits = [m for m in tar_members_by_range(read_range) if _is_chrom_file(m[0], chrom)]
    if not hits:
        raise FileNotFoundError(f"no chr{chrom} file in the archive")
    _, size, offset = hits[0]
    return slice_regenie_gz(read_range(offset, size), start38, end38, out_path)


def http_range_reader(url, retries=5):
    """Byte-range reader for a URL whose server honours Range requests (e.g. S3)."""
    def read_range(start, length):
        req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{start + length - 1}"})
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    if r.status != 206:
                        raise IOError(f"HTTP {r.status}: the server ignored the Range header")
                    return r.read()
            except (urllib.error.URLError, TimeoutError, ConnectionError):
                if attempt == retries - 1:
                    raise
                time.sleep(2 ** attempt)
    return read_range


def cis_out_name(protein_id, gene, multi_gene):
    """Output file name. Assays mapped to several genes (e.g. IL12A_IL12B) get one
    cis file per gene, because each gene has its own cis window."""
    stem = protein_id.replace(":", "_")
    return f"{stem}.{gene}.cis.regenie.gz" if multi_gene else f"{stem}.cis.regenie.gz"


def presigned_url(syn, entity_id):
    ent = syn.get(entity_id, downloadFile=False)
    body = json.dumps({"requestedFiles": [{"fileHandleId": ent.dataFileHandleId,
                                           "associateObjectId": entity_id,
                                           "associateObjectType": "FileEntity"}],
                       "includePreSignedURLs": True, "includeFileHandles": False})
    res = syn.restPOST("/fileHandle/batch", body=body, endpoint=syn.fileHandleEndpoint)
    return res["requestedFiles"][0]["preSignedURL"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--folder-id", required=True, help="Synapse ID of the UKB-PPP European discovery folder")
    ap.add_argument("--protein-map", required=True)
    ap.add_argument("--proteins", help="text file of UKBPPP_ProteinIDs or gene symbols (default: all)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--window-kb", type=float, default=500)
    ap.add_argument("--workers", type=int, default=8, help="proteins fetched in parallel")
    ap.add_argument("--full-download", action="store_true",
                    help="download whole tars instead of byte ranges (needs a few GB free per worker)")
    ap.add_argument("--tmpdir", default=None, help="where tars go with --full-download")
    a = ap.parse_args()

    import synapseclient  # imported here so the rest of the module works without it

    syn = synapseclient.Synapse(silent=True)
    syn.login(authToken=os.environ.get("SYNAPSE_AUTH_TOKEN"), silent=True)
    syn_lock = threading.Lock()

    pm = pd.read_csv(a.protein_map, sep="\t", dtype={"chr": str})
    pm["gene"] = pm["HGNC.symbol"]
    multi = set(pm.UKBPPP_ProteinID[pm.UKBPPP_ProteinID.duplicated()])
    if a.proteins:
        wanted = {l.strip() for l in open(a.proteins) if l.strip()}
        pm = pm[pm.UKBPPP_ProteinID.isin(wanted) | pm["HGNC.symbol"].isin(wanted)]
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tmp = Path(a.tmpdir or outdir / "_tmp")

    children = {c["name"]: c["id"] for c in syn.getChildren(a.folder_id, includeTypes=["file"])}
    w = int(a.window_kb * 1000)

    def fetch(r):
        pid = r.UKBPPP_ProteinID
        out = outdir / cis_out_name(pid, r.gene, pid in multi)
        if out.exists():
            return "done", pid, None
        stem = pid.replace(":", "_")
        match = [n for n in children if n.startswith(stem) and n.endswith(".tar")]
        if not match:
            return "skip", pid, f"no tar named {stem}*.tar in folder"
        chrom, start, end = str(r.chr), max(0, int(r.gene_start) - w), int(r.gene_end) + w
        if a.full_download:
            tmp.mkdir(exist_ok=True)
            with syn_lock:
                ent = syn.get(children[match[0]], downloadLocation=str(tmp))
            try:
                n = cis_slice_from_tar(ent.path, chrom, start, end, out)
            finally:
                Path(ent.path).unlink(missing_ok=True)
        else:
            with syn_lock:
                url = presigned_url(syn, children[match[0]])
            n = cis_slice_by_range(http_range_reader(url), chrom, start, end, out)
        return "ok", pid, f"{n} cis variants"

    counts = {"ok": 0, "done": 0, "skip": 0, "fail": 0}
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures = {pool.submit(fetch, r): f"{r.UKBPPP_ProteinID} ({r.gene})" for r in pm.itertuples()}
        for i, f in enumerate(as_completed(futures), 1):
            try:
                status, pid, msg = f.result()
            except Exception as e:  # keep going; a re-run retries failures
                status, pid, msg = "fail", futures[f], f"{type(e).__name__}: {e}"
            counts[status] += 1
            if status in ("ok", "skip", "fail"):
                print(f"[{status}] {i}/{len(futures)} {futures[f]}: {msg}", flush=True)
    print(f"finished: {counts}", flush=True)


if __name__ == "__main__":
    main()
