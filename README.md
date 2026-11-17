# Proteome-wide cis-MR and colocalisation scan for Parkinson's disease

A Nextflow pipeline that tests every plasma protein in the UK Biobank Pharma Proteomics Project (UKB-PPP, ~2,900 Olink proteins) as a potential cause of Parkinson's disease. For each protein it uses variants near the protein's own gene (cis-pQTLs) as genetic instruments, estimates the causal effect on PD risk with Mendelian randomisation (MR), and tests with colocalisation whether the protein and PD share the same causal variant. Proteins are FDR-corrected, sorted into evidence tiers and summarised in an HTML report.

This extends my MSc dissertation, which assessed a single protein ([GPNMB](https://github.com/barbavegeta/GPNMB-Parkinson-MR-Colocalisation)), into the scalable target-prioritisation workflow used by human-genetics teams in drug discovery.

It has been run on both large plasma pQTL resources: the UKB-PPP European discovery set (2,940 Olink assays, about 34,000 people) and the 3,211 SomaScan measurements of Sun et al. 2018 (INTERVAL, N = 3,301), which the eQTL Catalogue distributes openly. UKB-PPP's 1.6 TB of summary statistics are read with HTTP byte-range requests, so the whole proteome comes down as under 100 GB (see [Running on real data](#running-on-real-data)).

## Results: UKB-PPP Olink plasma proteome

![MR estimates for 2,016 testable Olink proteins, coloured by evidence tier](docs/img/volcano_ukbppp.png)

Of 2,940 Olink assays (2,957 rows: 15 assays measure a protein encoded by two or three genes and are tested once per gene), 2,016 had a genome-wide significant cis-pQTL to test, 848 had none, and 90 sit on chromosome X, which the PD GWAS does not cover. Eight pass FDR < 0.05:

| Tier | Protein | Locus | OR per SD (95% CI) | MR q | Coloc |
|---|---|---|---|---|---|
| 1: MR + colocalisation | CTF1 | 16p11.2 (KAT8/BCKDK) | 8.07 (4.14-15.75) | 9×10⁻⁷ | PP.H4 0.95 |
| 1: MR + colocalisation | STX4 | 16p11.2 (KAT8/BCKDK) | 0.09 (0.04-0.21) | 1×10⁻⁵ | PP.H4 0.89 |
| 2: suggestive colocalisation | LRRC37A2 | 17q21.31 (MAPT inversion) | 0.83 (0.77-0.90) | 0.006 | PP.H4 0.58 |
| 3: distinct causal variants | SNCA | SNCA | 0.10 (0.05-0.21) | 9×10⁻⁷ | PP.H3 1.00 |
| 3: distinct causal variants | PRSS53 | 16p11.2 (KAT8/BCKDK) | 0.90 (0.87-0.94) | 9×10⁻⁵ | PP.H3 0.99 |
| 3: distinct causal variants | PRSS8 | 16p11.2 (KAT8/BCKDK) | 0.46 (0.33-0.64) | 0.002 | PP.H3 1.00 |
| 3: distinct causal variants | EFNA1 | GBA1 | 1.23 (1.11-1.36) | 0.031 | PP.H3 1.00 |
| 3: distinct causal variants | HLA-DRA | HLA | 0.80 (0.72-0.90) | 0.042 | PP.H3 0.98 |

All eight sit in windows with a genome-wide significant PD signal, at known PD risk loci. Colocalisation is what separates them. For five, the protein and PD signals come from different variants, so the MR association rides on LD with the PD variant rather than on the protein. CTF1 and STX4 do colocalise, but they are evidence for the 16p11.2 locus rather than for either gene. Each rests on one modest pQTL (lead p = 1×10⁻¹¹ and 2×10⁻⁸), the implied effects are implausibly large and point in opposite directions in the same gene-dense region, and PP.H4 falls to 0.65 and 0.45 under a stricter prior (p12 = 10⁻⁶). LRRC37A2 lies in the 17q21.31 inversion, where long-range LD makes gene-level attribution unreliable. So the scan finds PD loci without being told where they are, but no single protein comes out as a confident causal candidate. That is a realistic outcome for a proteome-wide scan, and the tiers and flags are there to make it visible.

Full results: [`docs/results_ukbppp/scan_results.tsv`](docs/results_ukbppp/scan_results.tsv), per-method MR estimates and the [HTML report](docs/results_ukbppp/pd_proteome_mr_report.html). The scan took about 40 minutes on a 12-core laptop, after a four-hour download.

![Forest plot of the FDR-significant Olink proteins](docs/img/forest_ukbppp.png)

### Olink against SomaScan

The SomaScan hits (below) mostly do not carry over to Olink:

| Protein | SomaScan (Sun et al. 2018) | UKB-PPP Olink |
|---|---|---|
| FCGR2A | Tier 1: OR 1.06, PP.H4 0.96 | No MR evidence: OR 1.05, p = 0.07, PP.H4 0.11 |
| GPNMB | Tier 2-3: OR 1.50-1.56 | No MR evidence: OR 1.08, p = 0.53, PP.H3 1.00 |
| CTSB | Tier 2: OR 0.82, PP.H4 0.73 | Same direction, below FDR: OR 0.87, p = 0.03, PP.H4 0.73 |
| PTGR1 | Tier 2: OR 0.84, PP.H4 0.59 | Same direction, below FDR: OR 0.85, p = 8×10⁻⁴ (q = 0.12), PP.H4 0.49 |
| SNCA | Tier 3 | Tier 3 |
| BAG3 | Tier 3 | No MR evidence, PP.H3 1.00 |

For GPNMB this matches my dissertation: the Olink pQTL and the PD signal are distinct, while SomaScan gives a positive MR estimate. Disagreement between the platforms is common in pQTL studies. SomaScan aptamers and Olink antibodies can detect different forms of a protein, and protein-altering variants can create binding artefacts on either platform. CTSB and PTGR1 agree in direction and size on both platforms, which makes them the most consistent candidates from the two scans.

## Results: SomaScan plasma proteome (Sun et al. 2018)

![MR estimates for 767 testable proteins, coloured by evidence tier](docs/img/volcano_sun2018.png)

Of 3,211 SomaScan aptamers, 767 (713 genes) had a genome-wide significant cis-pQTL that could be tested. The rest had no instrument at p < 5×10⁻⁸ in 3,301 people. Nine aptamers pass FDR < 0.05:

| Tier | Protein | OR per SD (95% CI) | MR q | Coloc PP.H4 |
|---|---|---|---|---|
| 1: MR + colocalisation | FCGR2A | 1.06 (1.03-1.09) | 0.004 | 0.96 |
| 2: suggestive colocalisation | GPNMB (5080-131) | 1.50 (1.29-1.75) | 4×10⁻⁵ | 0.56 |
| 2: suggestive colocalisation | GPNMB (8240-207) | 1.53 (1.30-1.80) | 4×10⁻⁵ | 0.54 |
| 2: suggestive colocalisation | CTSB | 0.82 (0.74-0.91) | 0.012 | 0.73 |
| 2: suggestive colocalisation | PTGR1 | 0.84 (0.76-0.92) | 0.037 | 0.59 |
| 3: distinct causal variants | SNCA (8458-111, 8458-16) | 0.40-0.42 | 2×10⁻²² | 0.00 |
| 3: distinct causal variants | GPNMB (8289-8) | 1.56 (1.31-1.85) | 8×10⁻⁵ | 0.38 |
| 3: distinct causal variants | BAG3 | 1.62 (1.30-2.03) | 0.003 | 0.01 |

Unsupervised, the scan lands on known Parkinson's disease loci: FCGR2A, CTSB, GPNMB, SNCA and BAG3 are all at PD GWAS risk loci. Two of the three GPNMB aptamers reproduce my dissertation's single-instrument result (OR 1.32, N = 996) in a sample 3.3 times larger. SNCA and BAG3 show why the colocalisation tier matters. Their MR p-values are among the strongest in the scan, but the pQTL and PD signals come from different variants (PP.H3 ≈ 1), so the MR association most likely reflects LD with the true PD variant rather than the protein. Most estimates rest on one instrument: 255 of the 767 tested aptamers had more than one.

Full results: [`docs/results_sun2018/scan_results.tsv`](docs/results_sun2018/scan_results.tsv) (all 3,211 aptamers), per-method MR estimates and the [HTML report](docs/results_sun2018/pd_proteome_mr_report.html). The run took 43 minutes on a 12-core laptop.

![Forest plot of the FDR-significant proteins](docs/img/forest_sun2018.png)

## What it does

```
UKB-PPP pQTLs (per protein) ─┐
                             ├─ cis window (gene ± 500 kb, GRCh38) → QC (INFO, MAF)
PD GWAS (ieu-b-7, GRCh37) ───┘        │
                                      ├─ harmonise alleles on GRCh37 position (strand flips, palindromes)
                                      ├─ instruments: p < 5e-8, F ≥ 10, LD-clumped (r² < 0.001) on 1000G EUR
                                      ├─ MR: Wald ratio / IVW (fixed + random effects), MR-Egger, weighted median
                                      └─ coloc ABF (p1 = p2 = 1e-4, p12 = 1e-5) + prior sensitivity
                                                       │
            all proteins → Benjamini–Hochberg FDR → evidence tiers → HTML report
```

| Tier | Rule | Meaning |
|---|---|---|
| 1 | MR q < 0.05 and PP.H4 ≥ 0.8 | Protein and PD most likely share one causal variant |
| 2 | MR q < 0.05 and PP.H4 0.5–0.8 | Supportive, not conclusive |
| 3 | MR q < 0.05 and PP.H3 ≥ 0.5 | Different causal variants: MR probably reflects LD, not the protein |
| 4 | MR q < 0.05, coloc inconclusive | Usually an underpowered outcome signal |

Each protein runs as an independent Nextflow task, so the scan parallelises across a laptop's cores or an HPC/cloud cluster.

## Validation against R

The statistics are implemented in Python (NumPy/SciPy) and checked against the reference R packages. The test suite reproduces the R results from my dissertation to at least 8 significant figures:

| Check | R result (dissertation) | This pipeline |
|---|---|---|
| `coloc.abf`, SomaScan GPNMB vs PD (150 SNPs) | PP.H4 = 0.943955 | PP.H4 = 0.943955 |
| `coloc.abf`, UKB-PPP Olink GPNMB vs PD (4,018 SNPs) | PP.H3 = 0.999804 | PP.H3 = 0.999804 |
| coloc prior sensitivity (7 values of p12) | matches | matches |
| IVW fixed + multiplicative random effects, Cochran's Q (3 clumping settings) | e.g. OR 1.079, RE p = 0.473 | identical |
| Wald ratio, lead variant rs75801644 | β = −0.06717, SE 0.04832 | identical |

Run them with `python -m pytest -q` (17 tests). `validation/cross_check_with_R.R` repeats the coloc comparison with `coloc::coloc.abf` for any protein from a pipeline run.

## Test run

```bash
pip install -r requirements.txt
nextflow run main.nf -profile test
```

The test data combine the **real GPNMB locus** (UKB-PPP Olink pQTL and ieu-b-7 PD GWAS) with four **simulated proteins**, each built to land in a known tier. They were simulated from a genotype panel so LD, clumping and colocalisation behave realistically:

| Protein | Scenario | Result |
|---|---|---|
| SYN_SHARED | one causal pQTL that also drives PD | Tier 1 (OR 1.32, PP.H4 0.99) |
| SYN_MULTI | three independent causal pQTLs, all mediated by the protein | Tier 1, 3 instruments (OR 1.42, PP.H4 1.00) |
| SYN_DISTINCT | PD signal from a nearby variant in LD | Tier 3 (MR p = 0.014, PP.H3 1.00) |
| SYN_NULL | pQTL with no effect on PD | No MR evidence (PP.H1 0.99) |
| GPNMB (real) | see below | No MR evidence at the lead variant; flagged for multi-signal follow-up |

Outputs: `results_test/scan_results.tsv`, `results_test/pd_proteome_mr_report.html` and `pipeline_info/`. A copy is in [`docs/example_results/`](docs/example_results/).

![GPNMB regional plot](docs/img/gpnmb_regional.png)

GPNMB shows why the tiers and flags are there. The strongest Olink pQTL is a low-frequency variant (rs75801644) with no PD association, so lead-variant MR is null and single-variant colocalisation strongly prefers distinct signals (PP.H3 = 0.9998). But the pQTL has several independent signals, and PD's association peak lies over them. In a follow-up to my dissertation, SuSiE-based colocalisation found a secondary pQTL signal overlapping the PD signal, though which variant leads it depends on the LD reference panel. The pipeline flags such loci (`H3_at_PD_locus_check_multi_signal`) instead of discarding them.

## Running on real data

**SomaScan (Sun et al. 2018), no login needed.** Download the eQTL Catalogue file for dataset QTD000584 (`ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000035/QTD000584/QTD000584.all.tsv.gz`, 730 MB), then:

```bash
python bin/prepare_outcome.py --gwas ieu-b-7.vcf.gz --outdir outcome_by_chr
python scripts/prepare_eqtl_catalogue.py --sumstats QTD000584.all.tsv.gz \
    --outcome-dir outcome_by_chr --outdir data/sun2018_cis --samplesheet samplesheet_sun2018.csv
nextflow run main.nf --samplesheet samplesheet_sun2018.csv --exposure_format eqtlcat \
    --outcome ieu-b-7.vcf.gz --ld_ref 'ref/EUR_chr{chr}' --outdir results_sun2018
```

`prepare_eqtl_catalogue.py` works for any eQTL Catalogue dataset. It looks up GRCh37 positions through rsIDs in the outcome GWAS and fetches gene coordinates from the Ensembl REST API.

**UKB-PPP (Olink).**

1. **PD GWAS.** Download `ieu-b-7.vcf.gz` (Nalls et al. 2019, excluding 23andMe; 33,674 cases, 449,056 controls) from the IEU OpenGWAS project.
2. **LD reference.** 1000 Genomes Phase 3 EUR PLINK files (GRCh37), either as one genome-wide set or per chromosome (use `{chr}` in the prefix, e.g. `ref/EUR_chr{chr}`).
3. **pQTLs.** UKB-PPP European discovery summary statistics are on Synapse (folder `syn51365303`), one tar per protein, about 1.6 TB in total. `scripts/fetch_ukbppp_cis.py` never downloads a whole tar. The files sit in S3, which serves byte ranges, so it reads the tar headers (a few kB), downloads only the chromosome carrying the gene (~30 MB) and keeps the cis window. The whole proteome comes down as under 100 GB instead of 1.6 TB, several proteins in parallel (`--workers`), and a re-run skips proteins already done. A quick first pass is the proteins near PD risk loci:

```bash
python bin/prepare_outcome.py --gwas ieu-b-7.vcf.gz --outdir outcome_by_chr
python scripts/select_proteins_near_pd_loci.py --outcome-dir outcome_by_chr \
    --protein-map assets/olink_protein_map_3k_v1.tsv --out proteins_pd_loci.txt
export SYNAPSE_AUTH_TOKEN=...            # after accepting the UKB-PPP data terms
python scripts/fetch_ukbppp_cis.py --folder-id syn51365303 \
    --protein-map assets/olink_protein_map_3k_v1.tsv --proteins proteins_pd_loci.txt \
    --outdir data/ukbppp_cis --workers 16
python scripts/build_samplesheet.py --protein-map assets/olink_protein_map_3k_v1.tsv \
    --sumstats-dir data/ukbppp_cis --out samplesheet.csv

nextflow run main.nf --samplesheet samplesheet.csv --outcome ieu-b-7.vcf.gz \
    --ld_ref ref/EUR_chr{chr} --outdir results_pd
```

Leave out `--proteins` to fetch the whole proteome. Use `-profile slurm` on a cluster, `-profile docker` (build with `docker build -t pd-proteome-mr-scan:0.1.0 .`) or `-profile conda`. Every threshold can be set on the command line (`--window_kb`, `--p_threshold`, `--clump_r2`, `--p12`, `--fdr`; see `nextflow.config`).

## Repository layout

```
main.nf, nextflow.config     pipeline and profiles (test, docker, conda, slurm)
bin/                         CLI steps called by Nextflow
src/pdmr/                    library: io, harmonise, clump, mr, coloc, ranking, plots, report
scripts/                     data preparation (UKB-PPP fetch, eQTL Catalogue split, samplesheet, PD-locus selection, test data)
tests/                       validation against R, unit and end-to-end tests, fixtures
validation/                  R cross-check script
assets/                      UKB-PPP Olink protein map (gene coordinates, GRCh38)
docs/results_ukbppp/         real-data results on the UKB-PPP Olink proteome
docs/results_sun2018/        real-data results on the Sun et al. 2018 SomaScan proteome
```

## Limitations

- Colocalisation assumes one causal variant per trait in the region. Loci with several pQTL signals need SuSiE-based colocalisation. For UKB-PPP, SuSiE with a 503-person 1000 Genomes LD panel proved unstable (see the [GPNMB repository](https://github.com/barbavegeta/GPNMB-Parkinson-MR-Colocalisation)). For Sun et al. 2018, the eQTL Catalogue publishes SuSiE fine-mapping with in-sample LD, which is the planned route for multi-signal colocalisation.
- In gene-dense regions (16p11.2, 17q21.31, HLA) several proteins can colocalise with the same PD signal. Colocalisation then points to the locus, not to a gene.
- cis-MR assumes the instruments affect PD only through the protein. Protein-altering variants can change how the Olink antibodies bind (epitope effects) without changing protein level.
- SomaScan aptamers can bind more than their intended protein, and a pQTL can reflect altered aptamer binding rather than abundance. The 3,301-person sample also leaves most proteins without a cis instrument.
- The results are genetic prioritisation evidence, not proof of causality or druggability.

## Data and citations

UKB-PPP: Sun BB et al. *Nature* 2023. SomaScan: Sun BB et al. *Nature* 2018, via the eQTL Catalogue (Kerimov N et al. *Nat Genet* 2021). PD GWAS: Nalls MA et al. *Lancet Neurol* 2019 (ieu-b-7 via IEU OpenGWAS). Colocalisation: Giambartolomei C et al. *PLoS Genet* 2014. MR methods: Burgess S et al. 2013; Bowden J et al. 2015, 2016. The test fixtures under `tests/fixtures/` are derived from these public summary statistics.
