#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
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
}

#' Differently process dfs from rmats results depending on AS type
#' @param .df dataframe from rMATs
#' @param .start coord start
#' @param .end coord end
#' @param as_type alternative splicing type
#' @param sj_type type of splice junction
#' @param FDR is the FDR cutoff
#' @return parsed and filter rmats output
#'
process_RMATS <- function(df, .start, .end, as_type, sj_type) {
  .df <- .df %>%
    dplyr::select(
      chr, !!.start, !!.end, strand, comparison, FDR, IncLevelDifference
    ) %>%
    dplyr::rename(c(start = !!.start, end = !!.end))
  if (endsWith(as_type, "SS")) {
    .df <- .df %>%
      mutate(start = pmin(start, end), end = pmax(start, end))
  }
  .df$sj_type <- sj_type
  .df
}

files_mat <- str_split(files, "/", simplify = T)
comparison <- files_mat[, ncol(files_mat) - 1]
as_type <- gsub(".MATS.JC.txt", x = files_mat[, ncol(files_mat)], "")

dfs <- lapply(files, read.table, header = 1)
dfs <- lapply(
  seq_along(dfs),
  function(i) {
    cbind(
      dfs[[i]],
      comparison = comparison[i],
      as_type = as_type[i]
    )
  }
)
coordinate_mapping <- read.csv(
  text = "as_type,.start,.end,sj_type
SE,upstreamEE,downstreamES,ES_SJ
SE,upstreamEE,exonStart_0base,EI_SJ
IR,upstreamEE,downstreamES,IR
A5SS,longExonEnd,flankingES,A5SS_SSJ
A5SS,shortEE,flankingES,A5SS_ISJ
A3SS,flankingEE,longExonStart_0base,A3SS_SSJ
A3SS,flankingES,shortES,A3SS_ISJ"
)
# readability counts
res <- list()
for (.df in dfs) {
  .as_type <- unique(.df$as_type)
  message("Processing files from ", .as_type)
  .coordinate_mapping <- subset(
    coordinate_mapping, as_type == .as_type
  )
  for (i in rownames(.coordinate_mapping)) {
    .args <- as.list(.coordinate_mapping[i, ])
    .args$df <- .df

    res[[i]] <- base::do.call(process_RMATS, .args)
  }
}
# from https://stackoverflow.com/a/54075410/1694714
pipe_message <- function(.data, status) {
  message(status)
  .data
}

bind_rows(res) %>%
  pipe_message(str_glue("Number of SJ output by rMATs {nrow(.)}")) %>%
  filter(FDR < opt$cutoff) %>%
  pipe_message(
    str_glue(
      "Number of SJ output by rMATs after
    filtering (cutoff={opt$cutoff}) {nrow(.)}"
    )
  ) %>%
  mutate(chr = gsub(x = .$chr, replacement = "", pattern = "chr")) %>%
  arrange(FDR) %>%
  distinct(chr, start, end, strand, sj_type, .keep_all = TRUE) %>%
  pipe_message(
    str_glue(
      "Number of distinct SJ output by rMATs {nrow(.)}"
    )
  ) %>%
  write_csv(opt$out)