#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(dplyr)
  library(optparse)
  library(GenomicRanges)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "junctionseq/analysis/*_sigGenes.results.txt.gz",
    help = "Path with glob character to JunctioSeq result files.
    [default %default]", ,
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "junctionseq/junctionseq_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    type = "double",
    default = 0.05,
    help = "Discard junctions not called at the cutoff of [default %default]"
  )
)

if (exists("snakemake")) {
  opt <- list(
    input = snakemake@input,
    output = snakemake@output[[1]],
    cutoff = snakemake@params[["cutoff"]]
  )
  files <- opt$input
  file_names <- gsub(
    x = opt$input,
    pattern = "junctionseq/analysis/(.+)_sigGenes.results.txt.gz",
    replacement = "\\1"
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
  files <- Sys.glob(opt$input)
  file_names <- gsub(
    x = files,
    replacement = "\\1",
    pattern = sub(x = opt$input, pattern = "\\*", replacement = "(.*)")
  )
}


message("Loading JunctionSeq result")

read_junctionseq_out <- function(x) {
  tmp <- read.table(
    x,
    header = 1
  ) %>%
    filter(tmp$testable == T)
}

message("Loading processing the table")
res <- lapply(files, read_junctionseq_out)
names(res) <- file_names

res <- bind_rows(res, .id = "comparison")
message("Computing SJ ranks")

message(
  "Number of junctions output by JunctionSeq ",
  nrow(res)
)
res <- res %>%
  filter(padjust < opt$cutoff) %>%
  select(
    comparison,
    chr,
    start,
    end,
    strand,
    padjust,
    contains("log2FC"),
    geneID,
    featureType,
    expr_ref,
    expr_alt
  ) %>%
  mutate(method = "JunctionSeq")
message("Number of junctions after filtering ", nrow(res))
write_csv(res, opt$output)