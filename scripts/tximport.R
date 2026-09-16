#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[i + 1]
}

sample_string <- get_arg("--samples")
tx2gene_file <- get_arg("--tx2gene")
salmon_dir <- get_arg("--salmon-dir")
outdir <- get_arg("--outdir")

samples <- strsplit(sample_string, ",", fixed = TRUE)[[1]]

files <- file.path(salmon_dir, samples, "quant.sf")
names(files) <- samples

missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop(paste("Missing quant.sf:", paste(missing, collapse = ", ")))
}

tx2gene <- read.delim(
  tx2gene_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(tx2gene) < 2) {
  stop("tx2gene must contain at least two columns: transcript_id and gene_id")
}

tx2gene <- tx2gene[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

############################################################
# CONVERT MATRICES TO DATA FRAMES WITH EXPLICIT GENE COLUMN
############################################################

counts_df <- data.frame(
  GeneID_Ver = rownames(txi$counts),
  txi$counts,
  check.names = FALSE
)

tpm_df <- data.frame(
  GeneID_Ver = rownames(txi$abundance),
  txi$abundance,
  check.names = FALSE
)

length_df <- data.frame(
  GeneID_Ver = rownames(txi$length),
  txi$length,
  check.names = FALSE
)


############################################################
# WRITE OUTPUT
############################################################

write.table(
  counts_df,
  file = file.path(outdir, "gene_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tpm_df,
  file = file.path(outdir, "gene_tpm.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  length_df,
  file = file.path(outdir, "gene_length.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
