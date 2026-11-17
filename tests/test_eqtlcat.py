"""The eQTL Catalogue reader keeps the GRCh38 cis window and applies INFO / MAF QC."""
import pandas as pd

from pdmr.io import READERS, read_eqtlcat_slice


def test_eqtlcat_slice_window_and_qc(tmp_path):
    rows = pd.DataFrame({
        "chr": ["7", "7", "7", "7", "8"],
        "pos38": [1_000, 2_000, 3_000, 9_000, 2_000],
        "pos37": [1_100, 2_100, 3_100, 9_100, 2_100],
        "rsid": ["rs1", "rs2", "rs3", "rs4", "rs5"],
        "variant_id": ["a", "b", "c", "d", "e"],
        "ea": ["A", "C", "G", "T", "A"], "oa": ["G", "T", "A", "C", "G"],
        "eaf": [0.2, 0.005, 0.3, 0.4, 0.2],   # rs2 fails MAF
        "beta": [0.1, 0.2, 0.3, 0.4, 0.5], "se": [0.01] * 5,
        "p": [1e-3] * 5, "log10p": [3.0] * 5, "n": [3301] * 5,
        "info": [0.95, 0.95, 0.5, 0.95, 0.95],  # rs3 fails INFO
    })
    f = tmp_path / "trait.tsv"
    rows.to_csv(f, sep="\t", index=False)

    out = read_eqtlcat_slice(f, "7", 500, 5_000, info_min=0.8, maf_min=0.01)
    assert out.rsid.tolist() == ["rs1"]           # rs4 outside the window, rs5 on chr8
    assert out.pos37.tolist() == [1_100]
    assert READERS["eqtlcat"] is read_eqtlcat_slice
