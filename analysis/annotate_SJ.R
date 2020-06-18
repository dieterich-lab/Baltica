#!/usr/bin/env Rscript
# Title     : annotate_SJ from many DJU methods
# Objective : To annotate SJ called DS by many DJU methods with gene and transcript level information
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
                                 library(optparse)
                                 library(GenomicRanges)
                                 library(rtracklayer)
                               })

source("utils.R", chdir = TRUE)

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result", ,
    metavar = "character"
  ),
  make_option(
    c("-a", "--annotation"),
    type = "character",
    help = "Path to annotation (GTF/GFF format",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Path to output",
    metavar = "character"
  )
)
# enabling both cli or snakemake input
tryCatch({
           opt <- list(
             input = snakemake@input[[1]],
             annotation = snakemake@input[[2]],
             output = snakemake@output[[1]]
           )
         }, error = function() {
  opt <- parse_args(OptionParser(option_list = option_list))
})

message('Loading input')
data <- read.csv(opt$input)
gr <- GRanges(data)
gtf <- rtracklayer::import.gff2(opt$annotation)

message('Processing annotation')
tx <- subset(gtf, type == 'transcript')
tx <- subsetByOverlaps(tx, gr)
gtf <- gtf[gtf$transcript_id %in% unique(tx$transcript_id)]

ex_tx <- filter_multi_exon(gtf)
introns <- get_introns(ex_tx)

hits <- filter_hits_by_diff(gr, introns)

message('Annotating set of introns')
exon_number <- get_exon_number(ex_tx)
txid_to_gene <- setNames(tx$gene_name, nm = tx$transcript_id)
txid_to_tx_name <- setNames(tx$cmp_ref, nm = tx$transcript_id)
txid_to_classcode <- setNames(tx$class_code, nm = tx$transcript_id)

mcols(introns)$gene_name <- txid_to_gene[names(introns)]
mcols(introns)$transcript_name <- txid_to_tx_name[names(introns)]
mcols(introns)$class_code <- txid_to_classcode[names(introns)]
mcols(introns)$exon_number <- exon_number$value
introns$exon_number <- apply(introns$exon_number, 1, paste, collapse = '-')
names(introns) <- NULL

introns <- aggregate_metadata(introns)
hits <- filter_hits_by_diff(gr, introns, max_start = 2, max_end = 2)

gr_metadata <- data[queryHits(hits), c('comparison', 'method')]
gr_metadata$idx <- subjectHits(hits)
gr_metadata_agg <- aggregate(. ~ idx, gr_metadata, paste, collapse = ';')

message('Matching SJ to introns')

introns_with_match <- merge(
  as.data.frame(introns),
  gr_metadata_agg,
  by.x = "row.names",
  by.y = "idx"
)
introns_with_match <- subset(introns_with_match, select = -Row.names)
write.csv(introns_with_match, opt$output)
introns_wo_match  <- data[setdiff(seq_along(data), unique(queryHits(hits))) ]
write.csv(introns_wo_match, append_to_basename(opt$output))