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
    default = "rmats/*/*.MATS.JC.txt",
    help = "Path with glob character to Leafcutter result files.
    [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "rmats/rmats_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    type = "double",
    default = 0.05,
    help = "Discard junction with FDR < than --cutoff [default %default]"
  )
)
# enabling both baltica or snakemake input
if (exists("snakemake")) {
  opt <- list(
    input = snakemake@input,
    output = snakemake@output[[1]],
    cutoff = snakemake@params[["cutoff"]]
  )
  files <- opt$input
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
  files <- Sys.glob(opt$input)
  source("~/Baltica/scripts/utils.R")
}

#' Process for alternative splice site rMATs output files
#' @param df dataframe from rMATs
#' @param start name of the start column
#' @param end name of the end column
#' @param type flag for splice junction type
#' @param FDR is the FDR cutoff
#' @return GenomicRange of the selected SJ
#' @export
#'
process_RMATS_ass <- function(df, start = "flankingES", end = "shortES", type, FDR = 0.05) {
  df <- df %>%
    dplyr::filter(FDR < !!FDR) %>%
    dplyr::select(chr, !!start, !!end, strand, comparison, FDR, IncLevelDifference) %>%
    dplyr::rename(c(start = !!start, end = !!end)) %>%
    mutate(start = pmin(start, end), end = pmax(start, end))

  df <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
  df <- df[width(df) > 1, ]
  df <- unique(df)
  mcols(df)["type"] <- type

  suppressWarnings(try(seqlevelsStyle(df) <- "Ensembl"))


  df
}


#' Process for exon skipping and intron retention rMATs output files
#' @param df dataframe from rMATs
#' @param start name of the start column
#' @param end name of the end column
#' @param type flag for splice junction type
#' @param FDR is the FDR cutoff
#' @return GenomicRange of the selected SJ
#' @export
#'
process_RMATS <- function(df, start, end, type, FDR = 0.05) {
  df <- df %>%
    dplyr::filter(FDR < !!FDR) %>%
    dplyr::select(chr, !!start, !!end, strand, comparison, FDR, IncLevelDifference) %>%
    dplyr::rename(c(start = !!start, end = !!end))

  df <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = T)
  df <- unique(df)
  df$type <- type
  suppressWarnings(try(seqlevelsStyle(df) <- "Ensembl"))

  df
}

res <- split(files, str_split(files, "/", simplify = T)[, 2])

get_rmats_coord <- function(.files, .comparison) {
  message("Processing files for ", .comparison)
  x <- suppressWarnings(
    lapply(.files, readr::read_delim, "\t", col_types = c(
      .default = col_double(),
      GeneID = col_character(),
      geneSymbol = col_character(),
      chr = col_character(),
      strand = col_character(),
      IJC_SAMPLE_1 = col_character(),
      SJC_SAMPLE_1 = col_character(),
      IJC_SAMPLE_2 = col_character(),
      SJC_SAMPLE_2 = col_character(),
      IncLevel1 = col_character(),
      IncLevel2 = col_character()
    ))
  )
  names(x) <- str_split(.files, pattern = "[/.]", simplify = T)[, 3]
  for (i in names(x)) {
    x[[i]]$comparison <- .comparison
  }

  es_ssj <- process_RMATS(
    x$SE, "upstreamEE", "downstreamES", "ES_SJ",
    FDR = opt$cutoff
  )
  es_isj <- process_RMATS(
    x$SE, "upstreamEE", "exonStart_0base", "EI_SJ",
    FDR = opt$cutoff
  )
  ir <- process_RMATS(
    x$RI, "upstreamEE", "downstreamES", "IR",
    FDR = opt$cutoff
  )
  a5ss_SSJ <- process_RMATS_ass(
    x$A5SS, "longExonEnd", "flankingES", "A5SS_SSJ",
    FDR = opt$cutoff
  )
  a5ss_ISJ <- process_RMATS_ass(
    x$A5SS, "shortEE", "flankingES", "A5SS_ISJ",
    FDR = opt$cutoff
  )
  a3ss_SSJ <- process_RMATS_ass(
    x$A3SS, "flankingEE", "longExonStart_0base", "A5SS_SSJ",
    FDR = opt$cutoff
  )
  a3ss_ISJ <- process_RMATS_ass(
    x$A3SS, "flankingES", "shortES", "A5SS_ISJ",
    FDR = opt$cutoff
  )

  gr <- c(es_ssj, es_isj, ir, a5ss_SSJ, a5ss_ISJ, a3ss_SSJ, a3ss_ISJ)
  gr$method <- "rmats"

  gr
}

message("Loading and procesing rMATs files")
res <- lapply(setNames(names(res), names(res)), function(x) {
  get_rmats_coord(res[[x]], x)
})
res <- as(res, "GRangesList")
res <- unlist(res)

message("Number of junctions after filtering ", length(res))
write_csv(data.frame(res), opt$output)