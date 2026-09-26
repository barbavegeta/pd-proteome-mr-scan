# Cross-check one protein's Python results against the reference R packages.
#
# Usage (from the repository root, after a pipeline run):
#   Rscript validation/cross_check_with_R.R <work-dir-with-KEY.region.tsv.gz> <KEY> [exposure_N]
#
# Needs: install.packages("coloc"); data.table. Compares coloc::coloc.abf on the
# exact harmonised region the pipeline used with the PP.H0-H4 in KEY.summary.json.
# The golden tests in tests/test_validation_against_r.py already do this for the
# dissertation's GPNMB inputs; this script lets you repeat it for any hit.
suppressPackageStartupMessages({library(coloc); library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; key <- args[2]; n_exp <- if (length(args) > 2) as.numeric(args[3]) else 33000
n_case <- 33674; n_ctrl <- 449056

reg <- fread(file.path(dir, paste0(key, ".region.tsv.gz")))
maf <- pmin(reg$eaf_exp, 1 - reg$eaf_exp)
snp <- paste(reg$chr, reg$pos37, reg$ea_exp, reg$oa_exp, sep = ":")
res <- coloc.abf(
  dataset1 = list(beta = reg$beta_exp, varbeta = reg$se_exp^2, snp = snp, type = "quant", N = n_exp, MAF = maf),
  dataset2 = list(beta = reg$beta_out, varbeta = reg$se_out^2, snp = snp, type = "cc",
                  N = n_case + n_ctrl, s = n_case / (n_case + n_ctrl)),
  p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
py <- fromJSON(file.path(dir, paste0(key, ".summary.json")))
cmp <- data.frame(hypothesis = paste0("PP.H", 0:4),
                  R = as.numeric(res$summary[2:6]),
                  python = unlist(py[paste0("coloc_PP.H", 0:4)]))
cmp$abs_diff <- abs(cmp$R - cmp$python)
print(cmp, row.names = FALSE)
cat(if (max(cmp$abs_diff) < 1e-6) "\nMATCH\n" else "\nDIFFERENCE: check the exposure N passed as argument 3\n")
