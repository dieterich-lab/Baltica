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
    help = "Path with glob character to JunctioSeq result files. [default %default]",,
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "junctionseq/junctionseq_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

files <- Sys.glob(opt$input)

message("Loading JunctionSeq result")
add_containers <-  function (junctionseq_result){
    gr <- GRanges(junctionseq_result)
    hits <- findOverlaps(gr, drop.redundant=F, drop.self=F, ignore.strand=F)
    hits_groups <- split(gr[subjectHits(hits)], queryHits(hits))
    containers <- unlist(range(hits_groups))
    hits <- findOverlaps(gr, containers, type='within', select = 'all')
    hits.by.row <- split(gr[subjectHits(hits)], queryHits(hits))
    junctionseq_result$container <- as.character(unlist(range(hits.by.row)))
    junctionseq_result
}

read_junctionseq_out <- function(x) {
  read_table2(
    x,
    col_types = cols(
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
    )
  )
}

message("Loading processing the table")

res <- lapply(files, read_junctionseq_out)
names(res) <- gsub(
  x = files,
  replacement = '\\1',
  pattern = sub('\\*', '(.*)', files)
)

res  <- lapply(res, add_containers)
res <-  bind_rows(res, .id = 'comparison')

res <- res %>%
    arrange(comparison, container, expr_ref) %>%
    group_by(comparison, container) %>%
    mutate(ref_rank = rank(expr_ref, ties.method = "average")) %>%
    ungroup() %>%
    arrange(comparison, container, expr_alt) %>%
    group_by(comparison, container) %>%
    mutate(alt_rank = rank(expr_alt, ties.method = "average")) %>%
    ungroup()

res <- res %>% select(
    comparison,
    chr,
    start,
    end,
    strand,
    padjust,
    contains('log2FC('),
    geneID,
    featureType,
    container,
    expr_ref,
    expr_alt,
    ref_rank,
    alt_rank)

write_csv(res, opt$output)