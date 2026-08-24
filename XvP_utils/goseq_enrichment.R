#!/usr/bin/env Rscript
#
# GC-bias-corrected GO enrichment (goseq; Young et al. 2010, Genome Biology 11:R14).
#
# Deliberately minimal: this script does ONLY what requires R (goseq's probability
# weighting function and its Wallenius non-central hypergeometric test). Gene-set
# retrieval, universe construction, BH correction, term naming and classification all
# live on the Python side (XvP_utils/plotting/cross_comparison.py) so that every
# naming and threshold decision sits in one language.
#
# Usage:
#   Rscript goseq_enrichment.R <genes.csv> <gene2cat.csv> <out.csv> [<pwf.csv>]
#
#   genes.csv     header gene,de,bias  -- one row per gene in the universe;
#                 `gene` unique, `de` in {0,1}, `bias` numeric with no NA.
#   gene2cat.csv  header gene,category -- already subset to genes.csv$gene by Python.
#   out.csv       written here: category,p_wallenius,p_hypergeom,numDEInCat,numInCat
#   pwf.csv       optional; gene,bias,de,pwf for the P(DE | bias) diagnostic curve.
#
# Exit codes: 2 = bad arguments, 3 = goseq not installed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  cat("Usage: Rscript goseq_enrichment.R <genes.csv> <gene2cat.csv> <out.csv> [<pwf.csv>]\n",
      file = stderr())
  quit(status = 2L)
}
genes_csv <- args[1]; gene2cat_csv <- args[2]; out_csv <- args[3]
pwf_csv <- if (length(args) == 4L) args[4] else NA_character_

if (!requireNamespace("goseq", quietly = TRUE)) {
  cat("goseq is not installed for this R interpreter (", R.home(), ").\n", sep = "", file = stderr())
  quit(status = 3L)
}
suppressPackageStartupMessages(library(goseq))

g    <- read.csv(genes_csv,    stringsAsFactors = FALSE)
g2c  <- read.csv(gene2cat_csv, stringsAsFactors = FALSE)

# Named 0/1 vector is goseq's required input format; names must match gene2cat$gene.
de <- as.integer(g$de)
names(de) <- g$gene

# Probability weighting function: P(gene is DE | bias). This is the correction.
pwf <- nullp(de, bias.data = g$bias, plot.fit = FALSE)

# use_genes_without_cat = FALSE is goseq's default: the PWF is fit on the full
# universe, then genes belonging to no category are dropped from the term tests.
wal <- goseq(pwf, gene2cat = g2c)                              # GC-corrected
hyp <- goseq(pwf, gene2cat = g2c, method = "Hypergeometric")   # uncorrected reference

out <- merge(
  wal[, c("category", "over_represented_pvalue", "numDEInCat", "numInCat")],
  hyp[, c("category", "over_represented_pvalue")],
  by = "category", suffixes = c("_wal", "_hyp")
)
names(out)[names(out) == "over_represented_pvalue_wal"] <- "p_wallenius"
names(out)[names(out) == "over_represented_pvalue_hyp"] <- "p_hypergeom"
out <- out[, c("category", "p_wallenius", "p_hypergeom", "numDEInCat", "numInCat")]

write.csv(out, out_csv, row.names = FALSE)

if (!is.na(pwf_csv)) {
  pwf_out <- data.frame(gene = rownames(pwf), bias = pwf$bias.data,
                        de = pwf$DEgenes, pwf = pwf$pwf, stringsAsFactors = FALSE)
  write.csv(pwf_out, pwf_csv, row.names = FALSE)
}
