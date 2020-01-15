#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(tidylog)
})


option_list = list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "junctionseq/analysis/sigGenes.results.txt.gz",
    help = "Path to JunctionSeq result file",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "junctionseq/junctionseq_junctions.csv",
    help = "Output file",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# JunctionSeq outputs is one junction per line format
res <- read_table2(opt$input, col_types = cols(
  .default = col_double(),
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  featureID = col_character(),
  geneID = col_character(),
  countbinID = col_character(),
  testable = col_logical(),
  status = col_character(),
  allZero = col_logical(),
  strand = col_character(),
  transcripts = col_character(),
  featureType = col_character()
))

res <- res %>%
  filter(padjust < opt$pval_cuttoff, status == 'OK') %>%
  # filter_at(vars(contains('log2')), any_vars(abs(.) <= opt$foldchange_cutoff)) %>%
  select(chr, start, end, strand, padjust, contains('log2FC('), geneID, featureType) %>%
  pivot_longer(contains("log2FC("), names_to = "comparison", values_to = "l2fc") %>%
  filter(abs(l2fc) > opt$foldchange_cutoff)

comp_name <- res$comparison %>% str_match('log2FC\\((.*)\\)' )
res$comparison <- comp_name[,2 ] %>% str_replace('/', '_vs_')

write_csv(res, opt$output)


