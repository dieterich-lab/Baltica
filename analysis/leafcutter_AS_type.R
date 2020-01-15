#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(GenomicRanges)
  library(optparse)
})


option_list <- list(
  make_option(c("-i", "--input"), type="character", help=""),
  make_option(c("-a", "--annotation"), type="character", help=""),
  make_option(c("-o", "--output"), type="character", help="")
)


parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments = 1)
opt <- args$options
file <- args$args


x <- read_csv(
  opt$input,
  col_types = c(.default = col_number(),
                contrast = col_character(),
                chr = col_character(),
                intron = col_character(),
                strand = col_character(),
                cluster = col_character(),
                status = col_character(),
                genes = col_character()))

gtf <- rtracklayer::import.gff(opt$annotation)
exons <- subset(gtf, type == 'exon')
exons <- unique(exons)

gr <- GRanges(x)
gr.exons <- subsetByOverlaps(exons, reduce(gr))

as.type <- function(junction.gr, exon.gr) {
  if (as.logical(strand(junction.gr) != strand(exon.gr))) {
    warning('Junction and exon have different strand')
    return("NA")
  }

  if (length(junction.gr) != 1 | length(exon.gr) != 1) {
    warning('as.type function take one junction and exon per function call')
    return("NA")
  }

  # else if (j.strand == '+' & j.start %in% c(e.end, e.end + 1)) {
  j.strand <- as.character(strand(junction.gr)[1])
  j.start <- as.integer(start(junction.gr)[1])
  j.end <- as.integer(end(junction.gr)[1])
  e.start <- as.integer(start(exon.gr)[1])
  e.end <- as.integer(end(exon.gr)[1])

  is.annotated <- dplyr::case_when(
    j.strand == '+' & j.end == e.start ~ "JE",
    j.strand == '+' & j.start == e.end ~ "JS",
    j.strand == '-' & j.end == e.start ~ "JE",
    j.strand == '-' & j.start == e.end ~ "JS",
    TRUE ~ "NA"
  )

  if (is.annotated != "NA") {
    return(is.annotated)
  }

  if (j.strand == '+'){
    a = e.start - j.start
    b = j.end - e.end
    c = e.end - j.start
    d = j.end - e.start
  } else if (j.strand == '-'){
    a = j.end - e.end
    b = e.start - j.start
    c = j.end - e.start
    d = e.end - j.start
  } else {
    message('AS type definition needs strand information')
    return("NA")
  }

  type <- dplyr::case_when(
    all(c(a > 0, b > 0, c > 0, d > 0)) ~ "ES",
    all(c(a > 0, b < 0, c > 0, d > 0)) ~ "A5SS",
    all(c(a < 0, b < 0, c > 0, d > 0)) ~ "A3SS",
    TRUE ~ "NA")
  return(type)
}


j.exons <- apply(x, 1, function(row) {
  gr <- GRanges(row['junction'])
  subsetByOverlaps(exons, reduce(gr))
  }
)

junctions <- GRanges(x[['junction']])
split.junctions <- split(junctions, 1:length(junctions))

x$exons <- lapply(j.exons, function(x) { str_c(as.character(x), collapse=';') })

res <- c(); for (i in 1:length(junctions)){
  tmp <- lapply(
    as(j.exons[[i]], "GRangesList"),
    function(x){as.type(split.junctions[[i]], x)})

  tmp <- str_c(as.character(tmp), collapse=';')
  res <- c(res, tmp)
}
