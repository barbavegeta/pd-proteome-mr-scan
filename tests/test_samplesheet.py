"""build_samplesheet.py: assays encoded by several genes get one row and one cis file per gene."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from fetch_ukbppp_cis import cis_out_name  # noqa: E402


def test_multi_gene_assay_gets_one_row_per_gene(tmp_path):
    pm = pd.DataFrame({
        "UKBPPP_ProteinID": ["GPNMB:Q14956:OID20173:v1", "IL12A_IL12B:P29459_P29460:OID21327:v1",
                             "IL12A_IL12B:P29459_P29460:OID21327:v1"],
        "HGNC.symbol": ["GPNMB", "IL12A", "IL12B"],
        "chr": ["7", "3", "5"],
        "gene_start": [23246761, 159988835, 159314780],
        "gene_end": [23275108, 159996019, 159330487],
    })
    pm.to_csv(tmp_path / "map.tsv", sep="\t", index=False)
    d = tmp_path / "cis"
    d.mkdir()
    multi = set(pm.UKBPPP_ProteinID[pm.UKBPPP_ProteinID.duplicated()])
    names = [cis_out_name(p, g, p in multi) for p, g in zip(pm.UKBPPP_ProteinID, pm["HGNC.symbol"])]
    assert names == ["GPNMB_Q14956_OID20173_v1.cis.regenie.gz",
                     "IL12A_IL12B_P29459_P29460_OID21327_v1.IL12A.cis.regenie.gz",
                     "IL12A_IL12B_P29459_P29460_OID21327_v1.IL12B.cis.regenie.gz"]
    for n in names:
        (d / n).write_text("x")

    subprocess.run([sys.executable, str(ROOT / "scripts" / "build_samplesheet.py"), "--protein-map",
                    str(tmp_path / "map.tsv"), "--sumstats-dir", str(d), "--out", str(tmp_path / "s.csv")],
                   check=True, capture_output=True)
    s = pd.read_csv(tmp_path / "s.csv", dtype={"chr": str})
    assert s.protein_id.is_unique
    il12 = s[s.gene.str.startswith("IL12")].set_index("gene")
    assert il12.loc["IL12B", "protein_id"].endswith(":v1:IL12B")
    assert il12.loc["IL12B", "chr"] == "5"
    assert Path(il12.loc["IL12B", "exposure"]).name.endswith(".IL12B.cis.regenie.gz")
    assert Path(il12.loc["IL12A", "exposure"]).name.endswith(".IL12A.cis.regenie.gz")
