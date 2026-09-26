"""The tar-slicing and samplesheet helpers used when preparing real UKB-PPP data."""
import gzip
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from fetch_ukbppp_cis import cis_slice_from_tar  # noqa: E402

GPNMB = ROOT / "tests" / "fixtures" / "gpnmb_locus" / "gpnmb_ukbppp_chr7_22p8_23p8Mb_exposure.tsv"
PID = "GPNMB:Q14956:OID20173:v1"


def _fake_tar(tmp):
    d = tmp / "GPNMB_Q14956_OID20173_v1_Cardiometabolic"
    d.mkdir()
    for c in ("6", "7"):
        with open(GPNMB, "rb") as src, gzip.open(d / f"discovery_chr{c}_{PID}:Cardiometabolic.gz", "wb") as dst:
            shutil.copyfileobj(src, dst)
    tar = tmp / f"{d.name}.tar"
    with tarfile.open(tar, "w") as tf:
        tf.add(d, arcname=d.name)
    return tar


def test_cis_slice_and_samplesheet(tmp_path):
    tar = _fake_tar(tmp_path)
    cis_dir = tmp_path / "cis"
    cis_dir.mkdir()
    out = cis_dir / f"{PID.replace(':', '_')}.cis.regenie.gz"
    n = cis_slice_from_tar(tar, "7", 23_000_000, 23_500_000, out)
    d = pd.read_csv(out, sep=" ")
    assert n == len(d) > 0
    assert d.GENPOS.between(23_000_000, 23_500_000).all()

    sheet = tmp_path / "sheet.csv"
    subprocess.run([sys.executable, ROOT / "scripts" / "build_samplesheet.py",
                    "--protein-map", ROOT / "assets" / "olink_protein_map_3k_v1.tsv",
                    "--sumstats-dir", cis_dir, "--out", sheet], check=True)
    s = pd.read_csv(sheet)
    assert s.loc[0, "gene"] == "GPNMB" and str(s.loc[0, "chr"]) == "7"
    assert Path(s.loc[0, "exposure"]).exists()
