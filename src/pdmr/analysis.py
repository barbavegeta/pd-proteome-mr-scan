"""Per-protein cis-MR + colocalisation, from raw summary statistics to one result row."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .clump import LDReference, clump
from .coloc import Dataset, coloc_abf, prior_sensitivity
from .harmonise import harmonise
from .io import READERS, read_outcome_window
from .mr import run_mr

# MHC (GRCh38 chr6:28.51-33.48 Mb): long-range LD makes cis instruments unreliable.
MHC = ("6", 28_510_120, 33_480_577)

DEFAULTS = dict(window_kb=500, info_min=0.8, maf_min=0.01, p_threshold=5e-8, f_min=10.0,
                clump_r2=0.001, clump_kb=1000, palindrome_maf=0.42,
                p1=1e-4, p2=1e-4, p12=1e-5, outcome_ncase=33674, outcome_ncontrol=449056)


def analyse_protein(meta: dict, exposure_path, outcome_dir, ld_ref=None, exposure_format="ukbppp", **kw):
    """Run one protein. ``meta`` needs protein_id, gene, chr, start, end (GRCh38 gene bounds).

    Returns (summary dict, mr results frame, harmonised region frame, instruments frame).
    """
    cfg = {**DEFAULTS, **{k: v for k, v in kw.items() if v is not None}}
    chrom = str(meta["chr"])
    w = int(cfg["window_kb"] * 1000)
    start, end = int(meta["start"]) - w, int(meta["end"]) + w
    s = {"protein_id": meta["protein_id"], "gene": meta["gene"], "chr": chrom,
         "window_start38": max(0, start), "window_end38": end, "status": "ok", "flags": []}
    if chrom == MHC[0] and start < MHC[2] and end > MHC[1]:
        s["flags"].append("MHC")

    reader = READERS[exposure_format]
    exp = reader(exposure_path, chrom, max(0, start), end, info_min=cfg["info_min"], maf_min=cfg["maf_min"])
    s["n_cis_variants"] = int(len(exp))
    if exp.empty:
        s["status"] = "no_cis_variants"
        return s, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    out = read_outcome_window(outcome_dir, chrom, exp.pos37.min(), exp.pos37.max())
    h, audit = harmonise(exp, out, cfg["palindrome_maf"])
    s.update(audit)
    if h.empty:
        s["status"] = "no_overlap_with_outcome"
        return s, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    s["outcome_coverage"] = round(len(h) / len(exp), 3)

    lead_exp = exp.loc[exp.log10p.idxmax()] if "log10p" in exp else exp.loc[exp.p.idxmin()]
    s["exposure_lead_variant"] = lead_exp.get("variant_id", f"{chrom}:{lead_exp.pos37}")
    s["exposure_lead_log10p"] = float(lead_exp.log10p) if "log10p" in exp else float(-np.log10(lead_exp.p))
    if not ((h.pos37 == lead_exp.pos37) & (h.ea_exp == lead_exp.ea)).any():
        s["flags"].append("exposure_lead_missing_in_outcome")

    # ------------------------------------------------------------------ instruments
    h["F"] = (h.beta_exp / h.se_exp) ** 2
    cands = h[(h.p_exp < cfg["p_threshold"]) & (h.F >= cfg["f_min"])].copy()
    s["n_gw_significant"] = int(len(cands))
    inst = pd.DataFrame()
    s["instrument_selection"] = "none"
    if len(cands):
        ref = LDReference(ld_ref) if ld_ref else None
        if ref is not None and ref.available(chrom):
            c = cands.rename(columns={"ea_exp": "ea", "oa_exp": "oa", "p_exp": "p"})
            kept, n_missing = clump(c, ref, chrom, cfg["clump_r2"], cfg["clump_kb"])
            s["n_missing_from_ld_ref"] = n_missing
            if len(kept):
                inst = kept.rename(columns={"ea": "ea_exp", "oa": "oa_exp", "p": "p_exp"})
                s["instrument_selection"] = f"clumped r2<{cfg['clump_r2']}"
        if inst.empty:
            inst = cands.nlargest(1, "log10p") if "log10p" in cands else cands.nsmallest(1, "p_exp")
            s["instrument_selection"] = "lead variant (no LD reference match)"
    s["n_instruments"] = int(len(inst))

    mr_rows = []
    if len(inst):
        mr_in = inst.rename(columns={"beta_exp": "beta_exp", "se_exp": "se_exp",
                                     "beta_out": "beta_out", "se_out": "se_out"})
        primary, allres = run_mr(mr_in)
        mr_rows = allres
        s.update({f"mr_{k}": v for k, v in primary.items() if k in ("method", "beta", "se", "p", "OR", "OR_lci", "OR_uci")})
        if "Q_p" in primary:
            s["mr_Q_p"] = primary["Q_p"]
            if primary["Q_p"] < 0.05:
                s["flags"].append("heterogeneous_instruments")
        egger = next((r for r in allres if r["method"] == "MR-Egger"), None)
        if egger is not None:
            s["egger_intercept_p"] = egger["intercept_p"]
        if len(inst) > 1:
            s["flags"].append("multiple_independent_signals")
    else:
        s["status"] = "no_instrument"

    # ------------------------------------------------------------------ colocalisation
    cc = h[(h.se_exp > 0) & (h.se_out > 0) & h.eaf_exp.between(0, 1, inclusive="neither")]
    maf = np.minimum(cc.eaf_exp, 1 - cc.eaf_exp).values
    n_exp = float(np.nanmedian(cc.n_exp)) if "n_exp" in cc else float(np.nanmedian(cc.n))
    ntot = cfg["outcome_ncase"] + cfg["outcome_ncontrol"]
    d1 = Dataset(cc.beta_exp.values, cc.se_exp.values ** 2, "quant", N=n_exp, MAF=maf)
    d2 = Dataset(cc.beta_out.values, cc.se_out.values ** 2, "cc", N=ntot, s=cfg["outcome_ncase"] / ntot)
    co, snp_pp = coloc_abf(d1, d2, cfg["p1"], cfg["p2"], cfg["p12"], return_snp=True)
    s.update({f"coloc_{k}": v for k, v in co.items()})
    sens = prior_sensitivity(d1, d2, (1e-6, 1e-5, 1e-4), cfg["p1"], cfg["p2"])
    s["coloc_PP.H4_p12_1e-6"] = sens[1e-6]
    s["coloc_PP.H4_p12_1e-4"] = sens[1e-4]
    # H3 alongside a genome-wide significant PD signal is where single-variant coloc is least
    # reliable (e.g. a strong secondary pQTL may be the shared one): flag for SuSiE-coloc follow-up.
    s["outcome_min_p_in_window"] = float(cc.p_out.min())
    if co["PP.H3"] > 0.5 and cc.p_out.min() < 5e-8:
        s["flags"].append("H3_at_PD_locus_check_multi_signal")
    cc = cc.assign(snp_pp_h4=snp_pp)
    top = cc.loc[cc.snp_pp_h4.idxmax()]
    s["coloc_top_variant"] = top.get("rsid", f"{chrom}:{top.pos37}")
    s["coloc_top_snp_pp"] = float(top.snp_pp_h4)

    region_cols = [c for c in ["chr", "pos37", "pos38", "rsid", "variant_id", "ea_exp", "oa_exp", "eaf_exp",
                               "beta_exp", "se_exp", "p_exp", "beta_out", "se_out", "p_out"] if c in cc]
    region = cc[region_cols + ["snp_pp_h4"]]
    inst_cols = [c for c in ["chr", "pos37", "rsid", "variant_id", "ea_exp", "oa_exp", "beta_exp", "se_exp",
                             "p_exp", "F", "beta_out", "se_out", "p_out"] if c in inst]
    mr_df = pd.DataFrame(mr_rows)
    if len(mr_df):
        mr_df.insert(0, "gene", meta["gene"])
        mr_df.insert(0, "protein_id", meta["protein_id"])
    s["flags"] = ";".join(s["flags"])
    return s, mr_df, region, inst[inst_cols] if len(inst) else inst
