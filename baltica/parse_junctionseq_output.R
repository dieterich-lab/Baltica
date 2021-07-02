#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
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
  col_names <- "featureID
geneID
countbinID
testable
status
allZero
baseMean
baseVar
dispBeforeSharing
dispFitted
dispersion
pvalue
padjust
chr
start
end
strand
transcripts
featureType
padjust_noFilter
log2FC_alt_vs_ref
log2FCvst_alt_vs_ref
expr_ref
expr_alt
geneWisePadj"

  tmp <- read_table2(
    x,
    col_names = strsplit(col_names, "\n")[[1]],
    skip = 1,
    col_types = cols(
      .default = col_double(),
      chr = col_character(),
      featureID = col_character(),
      geneID = col_character(),
      countbinID = col_character(),
      testable = col_logical(),
      status = col_character(),
      allZero = col_logical(),
      strand = col_character(),
      transcripts = col_character(),
      featureType = col_character()
    )
  )

  tmp %>% filter(tmp$testable == T)
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