#!/usr/bin/env python3
"""Build the small test dataset used by ``-profile test`` and the pytest suite.

It combines:
* the real GPNMB cis region (UKB-PPP Olink pQTL + ieu-b-7 PD GWAS, chr7), and
* four simulated proteins on chr22 with a matching simulated LD panel, each
  built to exercise one scenario:
    SYN_SHARED    one causal pQTL that also drives PD      -> MR + coloc (H4)
    SYN_DISTINCT  PD signal from a nearby variant in LD     -> MR hit, coloc H3
    SYN_NULL      pQTL with no PD effect                    -> null MR, coloc H1
    SYN_MULTI     three independent causal pQTLs, mediated  -> multi-instrument IVW
Summary statistics are simulated from the panel's LD so clumping and
colocalisation behave realistically.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import to_bed
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "data"
REAL = ROOT / "tests" / "fixtures" / "gpnmb_locus"
rng = np.random.default_rng(2026)

N_IND, N_EXP = 503, 33000  # panel size (like 1000G EUR); pQTL sample size (like UKB-PPP)
N_CASE, N_CTRL = 33674, 449056
NEFF = (N_CASE + N_CTRL) * (N_CASE / (N_CASE + N_CTRL)) * (N_CTRL / (N_CASE + N_CTRL))
PAIRS = [("A", "G"), ("C", "T"), ("A", "C"), ("G", "T"), ("A", "T"), ("C", "G")]
PAIR_P = [0.3, 0.3, 0.17, 0.17, 0.03, 0.03]  # a few palindromes to exercise harmonisation

LOCI = {  # gene: (GRCh37 centre on chr22, scenario)
    "SYN_SHARED": (20_000_000, "shared"),
    "SYN_DISTINCT": (30_000_000, "distinct"),
    "SYN_NULL": (40_000_000, "null"),
    "SYN_MULTI": (45_000_000, "multi"),
}
SNPS_PER_LOCUS, BLOCK = 240, 30
OFFSET38 = -40_000  # arbitrary GRCh37 -> "GRCh38" shift so the two coordinate systems differ


def simulate_haplotype_block(n_hap, n_snp):
    """Latent AR(1) Gaussian per haplotype, thresholded -> alleles in LD that decays with distance."""
    rho = rng.uniform(0.85, 0.97)
    z = np.empty((n_hap, n_snp))
    z[:, 0] = rng.standard_normal(n_hap)
    for j in range(1, n_snp):
        z[:, j] = rho * z[:, j - 1] + np.sqrt(1 - rho ** 2) * rng.standard_normal(n_hap)
    freq = rng.uniform(0.05, 0.5, n_snp)
    thr = np.array([np.quantile(z[:, j], 1 - freq[j]) for j in range(n_snp)])
    return (z > thr).astype(np.int8)


def simulate_locus():
    blocks = [simulate_haplotype_block(2 * N_IND, BLOCK) for _ in range(SNPS_PER_LOCUS // BLOCK)]
    hap = np.hstack(blocks)
    g = hap[0::2] + hap[1::2]  # diploid dosage of the "effect" allele
    return g.astype(float)


def mvn(R, mean):
    L = np.linalg.cholesky(R + 1e-6 * np.eye(len(R)))
    return mean + L @ rng.standard_normal(len(R))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "exposures").mkdir(exist_ok=True)
    bim_rows, geno, vcf_rows, map_rows, sheet = [], [], [], [], []

    for gene, (centre, scen) in LOCI.items():
        g = simulate_locus()
        keep = g.std(0) > 0
        g = g[:, keep]
        m = g.shape[1]
        pos37 = np.sort(rng.choice(np.arange(centre - 400_000, centre + 400_000), m, replace=False))
        f = g.mean(0) / 2
        R = np.corrcoef(g.T)
        pairs = [PAIRS[i] for i in rng.choice(len(PAIRS), m, p=PAIR_P)]
        a0 = np.array([p[0] for p in pairs]); a1 = np.array([p[1] for p in pairs])
        centre_idx = np.argmin(np.abs(pos37 - centre))
        common = np.where((f > 0.15) & (f < 0.85))[0]

        b_exp = np.zeros(m); b_out = np.zeros(m)
        theta = 0.35  # log OR for PD per SD of protein in the mediated scenarios
        c = common[np.argmin(np.abs(common - centre_idx))]
        if scen == "shared":
            b_exp[c] = 0.09; b_out = theta * b_exp
        elif scen == "distinct":
            b_exp[c] = 0.09
            r = R[c]
            cand = [j for j in common if 0.3 < abs(r[j]) < 0.55]
            d = cand[0] if cand else common[(list(common).index(c) + 5) % len(common)]
            b_out[d] = 0.045 * np.sign(r[d])
        elif scen == "null":
            b_exp[c] = 0.09
        elif scen == "multi":
            picks = [c]
            for j in common[np.argsort(np.abs(common - centre_idx))]:
                if len(picks) == 3:
                    break
                if all(abs(R[j, k]) < 0.02 for k in picks):
                    picks.append(j)
            for k, eff in zip(picks, (0.08, 0.06, 0.05)):
                b_exp[k] = eff
            b_out = theta * b_exp

        z_exp = mvn(R, np.sqrt(N_EXP) * R @ b_exp)
        z_out = mvn(R, np.sqrt(NEFF) * R @ b_out)
        het = 2 * f * (1 - f)
        se_exp = 1 / np.sqrt(het * N_EXP)
        se_out = 1 / np.sqrt(het * NEFF)
        beta_exp, beta_out = z_exp * se_exp, z_out * se_out
        log10p_exp = -np.log10(2 * norm.sf(np.abs(z_exp)))
        lp_out = -np.log10(2 * norm.sf(np.abs(z_out)))

        # UKB-PPP style REGENIE file (BETA / A1FREQ refer to ALLELE1)
        exp = pd.DataFrame({
            "CHROM": 22, "GENPOS": pos37 + OFFSET38,
            "ID": [f"22:{p}:{x}:{y}:imp:v1" for p, x, y in zip(pos37, a0, a1)],
            "ALLELE0": a0, "ALLELE1": a1, "A1FREQ": f.round(6), "INFO": rng.uniform(0.85, 1, m).round(4),
            "N": N_EXP, "TEST": "ADD", "BETA": beta_exp, "SE": se_exp, "CHISQ": z_exp ** 2,
            "LOG10P": log10p_exp, "EXTRA": "NA"})
        path = OUT / "exposures" / f"discovery_chr22_{gene}.regenie.gz"
        exp.to_csv(path, sep=" ", index=False, compression="gzip")
        pid = f"{gene}:SIM000:OIDSIM:v1"
        s38 = int(pos37[centre_idx] + OFFSET38 - 5_000)
        map_rows.append({"UKBPPP_ProteinID": pid, "HGNC.symbol": gene, "chr": 22,
                         "gene_start": s38, "gene_end": s38 + 10_000, "ensembl_id": "", "Panel": "Simulated"})
        sheet.append({"protein_id": pid, "gene": gene, "chr": 22, "start": s38, "end": s38 + 10_000,
                      "exposure": f"exposures/{path.name}"})

        # Outcome VCF rows (ieu-b-7 style); a third of variants stored with alleles swapped
        swap = rng.random(m) < 0.33
        for i in range(m):
            ref, alt, es, af = (a1[i], a0[i], -beta_out[i], 1 - f[i]) if swap[i] else (a0[i], a1[i], beta_out[i], f[i])
            vcf_rows.append(f"22\t{pos37[i]}\trs9{gene[-3:].lower()}{i}\t{ref}\t{alt}\t.\tPASS\tAF={af:.4f}\t"
                            f"ES:SE:LP:AF:SS:NC:ID\t{es:.5f}:{se_out[i]:.5f}:{lp_out[i]:.5f}:{af:.4f}:"
                            f"{N_CASE + N_CTRL}:{N_CASE}:rs9{gene[-3:].lower()}{i}")
        for i in range(m):
            bim_rows.append(f"22\t22_{pos37[i]}\t0\t{pos37[i]}\t{a1[i]}\t{a0[i]}")
        geno.append(g)

    # LD panel (PLINK bed/bim/fam), chr22 only
    G = np.hstack(geno).astype("float32")
    ld = OUT / "ld_ref"
    ld.mkdir(exist_ok=True)
    to_bed(ld / "synthetic_eur_chr22.bed", G,
           properties={"chromosome": ["22"] * G.shape[1],
                       "sid": [r.split("\t")[1] for r in bim_rows],
                       "bp_position": [int(r.split("\t")[3]) for r in bim_rows],
                       "allele_1": [r.split("\t")[4] for r in bim_rows],
                       "allele_2": [r.split("\t")[5] for r in bim_rows]})

    # Outcome VCF: real GPNMB chr7 region + simulated chr22
    real = (REAL / "ieu_b_7_chr7_22p8_23p8Mb_PD_outcome.vcf").read_text().splitlines()
    with gzip.open(OUT / "pd_outcome_test.vcf.gz", "wt") as fh:
        for line in real:
            fh.write(line + "\n")
        for line in vcf_rows:
            fh.write(line + "\n")

    # Protein map and samplesheet (real GPNMB row first)
    gp = {"UKBPPP_ProteinID": "GPNMB:Q14956:OID20173:v1", "HGNC.symbol": "GPNMB", "chr": 7,
          "gene_start": 23235967, "gene_end": 23275108, "ensembl_id": "ENSG00000136235", "Panel": "Cardiometabolic"}
    pd.DataFrame([gp] + map_rows).to_csv(OUT / "protein_map_test.tsv", sep="\t", index=False)
    pd.DataFrame([{"protein_id": gp["UKBPPP_ProteinID"], "gene": "GPNMB", "chr": 7,
                   "start": gp["gene_start"], "end": gp["gene_end"],
                   "exposure": "../fixtures/gpnmb_locus/gpnmb_ukbppp_chr7_22p8_23p8Mb_exposure.tsv"}] + sheet
                 ).to_csv(OUT / "samplesheet_test.csv", index=False)
    print("test data written to", OUT)


if __name__ == "__main__":
    main()
