#!/usr/bin/env Rscript
# Title     : annotate_SJ from many DJU methods
# Objective : To annotate SJ called DS by many DJU methods with gene and
# transcript level information
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(readr)
  library(yaml)
  library(igraph)
  library(tidyverse)
})
#' Compute and filter hits based on the difference in the genomic start and end
#'
#' @param query first set of range
#' @param subject second set of ranges
#' @param max_start max_end absolute max difference at the start and end
#' coordinates, respectively
#' @return overlapping ranges given the contrain
#' @export
filter_hits_by_diff <- function(query, subject, max_start = 2, max_end = 2) {
  stopifnot(is(query, "GRanges"))
  stopifnot(is(subject, "GRanges"))
  hits <- findOverlaps(query, subject)
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  start_dif <- abs(start(query) - start(subject))
  end_dif <- abs(end(query) - end(subject))
  hits <- hits[start_dif <= max_start & end_dif <= max_end]
  hits
}

#' Creates exon by transcript list (ex_by_tx) and remove items
#' with a single exon
#'
#' @param gtf GRange file loaded with rtracklayer::import
#' @return  a list of exon by transcripts, excluding single exon trascripts
#' @export
filter_multi_exon <- function(gtf) {
  stopifnot(is(gtf, "GRanges"))

  ex <- subset(gtf, type == "exon")
  stopifnot(length(ex) > 0)

  # discard single exons transcripts
  multi_ex <- table(ex$transcript_id) > 1
  ex <- subset(ex, mcols(ex)$transcript_id %in% names(multi_ex[multi_ex]))
  ex_tx <- split(ex, ex$transcript_id)
  ex_tx
}

#' Compute a set of introns from GTF files
#'
#' @param gtf_path path to to the gtf file
#' @param read_gtf function to read the gtf_file
#' @return a GRange that obj with the introns named by their parent trancripts
#' @export
get_introns <- function(ex_tx) {
  stopifnot(is(ex_tx, "List"))

  introns <- psetdiff(unlist(range(ex_tx), use.names = FALSE), ex_tx)
  introns <- unlist(introns)
  introns
}

aggregate_annotation <- function(gr) {
  stopifnot(is(gr, "GRanges"))

  equal_hits <- findOverlaps(gr, type = "equal")
  equal_hits <- as.data.frame(equal_hits)
  equal_hits <- igraph::graph_from_data_frame(equal_hits)
  equal_hits_groups <- stack(igraph::groups(igraph::clusters(equal_hits)))


  gr_ <- as.data.frame(mcols(gr))
  gr_$coordinates <- as.character(gr)
  gr_[as.numeric(equal_hits_groups$values), "group"] <- equal_hits_groups$ind
  gr_ <- aggregate(. ~ group, gr_, unique)
  gr <- GenomicRanges::GRanges(gr_$coordinates)
  mcols(gr) <- gr_
  gr
}

#' Fetch exons pairs for an intron
#'
#' @param introns_by_transcript list of introns by transcript
#' @return a data.frame with acceptor and donor exon number for an intron
#' @export
get_exon_number <- function(ex_tx) {
  stopifnot(is(ex_tx, "List"))

  exon_n_by_transcript <- lapply(ex_tx, function(x) embed(x$exon_number, 2))
  exon_number <- tidyr::unnest(
    tibble::enframe(exon_n_by_transcript),
    cols = "value"
  )

  exon_number
}

default_input <- str_glue(
  "{method}/{method}_junctions.csv",
  method = c("junctionseq", "leafcutter", "majiq", "rmats")
) %>% str_c(collapse = ",")

#  ,/leafcutter_junctions.csv,majiq/majiq_junctions.csv,rmats/.csv"

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result",
    metavar = "character",
    default = default_input
  ),
  make_option(
    c("-a", "--annotation"),
    type = "character",
    help = "Path to annotation (GTF/GFF format)",
    metavar = "character",
    default = "stringtie/merged/merged.combined.gtf"
  ),
  make_option(
    c("-r", "--reference"),
    type = "character",
    help = "Path to reference annotation (GTF/GFF format)",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Path to output",
    metavar = "character",
    default = "results/SJ_annotated.csv"
  )
)
# enabling both baltica or snakemake input
if (exists("snakemake")) {
  opt <- list(
    input = paste0(
      c(
        snakemake@input[[1]],
        snakemake@input[[2]],
        snakemake@input[[3]],
        snakemake@input[[4]]
      ),
      collapse = ","
    ),
    reference = snakemake@params$ref,
    annotation = snakemake@input[[5]],
    output = snakemake@output[[1]]
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
}

files <- strsplit(opt$input, ",")[[1]]

if (!all(as.logical(lapply(files, file.exists)))) {
  stop("Input file not found.", call. = FALSE)
} else if (!file.exists(opt$annotation)) {
  stop("Annotation not found.", call. = FALSE)
}

gtf <- rtracklayer::import.gff2(opt$annotation)
.read_csv <- function(x) {
  suppressWarnings(
    readr::read_csv(x, col_types = readr::cols(
      chr = readr::col_character(),
      seqnames = readr::col_character()
    ))
  )
}
majiq_idx <- grep("majiq", files)
leafcutter_idx <- grep("leafcutter", files)
junctionseq_idx <- grep("junctionseq", files)
rmats_idx <- grep("rmats", files)

message("Processing DJU method output")
df <- list(
  majiq = .read_csv(files[[majiq_idx]]),
  leafcutter = .read_csv(files[[leafcutter_idx]]),
  junctionseq = .read_csv(files[[junctionseq_idx]]),
  rmats = .read_csv(files[[rmats_idx]])
)

gr <- suppressWarnings(c(
  GRanges(df$majiq),
  GRanges(df$leafcutter),
  GRanges(df$junctionseq),
  GRanges(df$rmats)
))

mcols(gr) <- NULL
mcols(gr)$method <- c(
  df$majiq$method,
  df$leafcutter$method,
  df$junctionseq$method,
  df$rmats$method
)
gr$method <- tolower(gr$method)
mcols(gr)$comparison <- c(
  df$majiq$comparison,
  df$leafcutter$comparison,
  df$junctionseq$comparison,
  df$rmats$comparison
)
mcols(gr)$score <- c(
  1 - df$majiq$probability_non_changing,
  1 - df$leafcutter$p.adjust,
  1 - df$junctionseq$padjust,
  1 - df$rmats$FDR
)

message("Processing de novo annotation")
tx <- subset(gtf, type == "transcript")
tx <- subsetByOverlaps(tx, gr)
gtf <- gtf[gtf$transcript_id %in% unique(tx$transcript_id)]

ex_tx <- filter_multi_exon(gtf)
introns <- get_introns(ex_tx)

if (!is.null(opt$reference)) {
  message("Processing reference annotation")

  reference <- rtracklayer::import.gff(opt$reference)
  ref_ex_tx <- filter_multi_exon(reference)

  ref_introns <- get_introns(ref_ex_tx)
  introns$is_novel <- !(introns %in% ref_introns)
}
message("preparing annotation")
exon_number <- get_exon_number(ex_tx)
txid_to_gene <- setNames(tx$gene_name, nm = tx$transcript_id)
txid_to_tx_name <- setNames(tx$cmp_ref, nm = tx$transcript_id)
txid_to_classcode <- setNames(tx$class_code, nm = tx$transcript_id)

mcols(introns)$gene_name <- txid_to_gene[names(introns)]
mcols(introns)$transcript_name <- txid_to_tx_name[names(introns)]
mcols(introns)$class_code <- txid_to_classcode[names(introns)]
mcols(introns)$exon_number <- exon_number$value
introns$exon_number <- apply(introns$exon_number, 1, paste, collapse = "-")

introns <- aggregate_annotation(introns)
introns$group <- NULL
introns_metadata <- mcols(introns)
introns_metadata <- introns_metadata %>%
  as_tibble() %>%
  rownames_to_column("metadata_group")

mcols(introns) <- NULL
names(introns) <- NULL
mcols(introns)$metadata_group <- seq_along(introns)

gr <- c(gr, introns)

message("Integration introns")

hits <- findOverlaps(gr, drop.redundant = T)
query <- gr[queryHits(hits)]
subject <- gr[subjectHits(hits)]
start_dif <- abs(start(query) - start(subject))
end_dif <- abs(end(query) - end(subject))
max_start <- 2
max_end <- 2
hits <- hits[start_dif <= max_start & end_dif <= max_end]

df_base <- as_tibble(mcols(gr)[to(hits), ])
df_base$group <- from(hits)

df <- df_base %>%
  group_by(
    group
  ) %>%
  fill(
    metadata_group,
    .direction = "up"
  ) %>%
  drop_na(method) %>%
  pivot_wider(
    id_cols = c(group, metadata_group),
    names_from = c(method, comparison),
    values_from = c(score),
    # if there is more than one intron entry for
    # the same comparison and method
    # use the highest scoring one
    # this happens for methods that test for multiple
    # AS types
    values_fn = c(score = max)
  )

message("Integrating annotation")
index <- split(gr[to(hits), ], from(hits))
group_to_annotation <- mcols(stack(index, "group"))
group_to_annotation <- group_to_annotation[c("group", "metadata_group")]
group_to_annotation <- group_to_annotation[
  !is.na(group_to_annotation$metadata_group),
]

index_regions <- stack(range(range(index)), "group")
index_regions$coord <- as.character(index_regions)
index_regions <- mcols(index_regions)
df <- plyr::join(df, as_tibble(group_to_annotation),
  match = "first",
  by = "group"
)
df <- plyr::join(df, introns_metadata, match = "first")
df <- plyr::join(df, as_tibble(index_regions), match = "first")
df <- df %>%
  mutate(coordinates = coalesce(coord)) %>%
  select(-coord, -group, -metadata_group) %>%
  mutate_if(
    is.list,
    list(~ map_chr(., . %>% str_c(collapse = ";")))
  )

readr::write_csv(df, opt$output)