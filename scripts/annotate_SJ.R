#!/usr/bin/env Rscript
# Title     : annotate_SJ from many DJU methods
# Objective : To annotate SJ called DS by many DJU methods with gene and transcript level information
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
  library(readr)
})

# from https://stackoverflow.com/a/15373917/1694714
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result", ,
    metavar = "character",
    default = "junctionseq/junctionseq_junctions.csv,leafcutter/leafcutter_junctions.csv,majiq/majiq_junctions.csv"
  ),
  make_option(
    c("-a", "--annotation"),
    type = "character",
    help = "Path to annotation (GTF/GFF format)",
    metavar = "character",
    default = "stringtie/merged/merged.combined.gtf"
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
if(exists('snakemake')){
  opt <- list(
    input = paste0(c(snakemake@input[[1]], snakemake@input[[2]], snakemake@input[[3]]), collapse = ','),
    annotation = snakemake@input[[4]],
    output = snakemake@output[[1]]
  )
  snakemake@source("utils.R")
  } else {
  opt <- parse_args(OptionParser(option_list = option_list))
  source(file.path(dirname(thisFile()), "utils.R"))
}
files <- strsplit(opt$input, ',')[[1]]
majiq_idx <- grep('majiq', files)
leafcutter_idx <- grep('leafcutter', files)
junctionseq_idx <- grep('junctionseq', files)

suppress_read_csv <- function(x) { suppressMessages(read_csv(x)) }

df  <-  list(
    majiq = suppress_read_csv( files[[majiq_idx]] ),
    leafcutter = suppress_read_csv( files[[leafcutter_idx]] ) ,
    junctionseq = suppress_read_csv( files[[junctionseq_idx]]  )
)

gr  <- c(
    GRanges(df$majiq),
    GRanges(df$leafcutter),
    GRanges(df$junctionseq)
)

mcols(gr) <- NULL
mcols(gr)$method <- c(df$majiq$method, df$leafcutter$method, df$junctionseq$method)
mcols(gr)$comparison <- c(df$majiq$comparison, df$leafcutter$comparison, df$junctionseq$comparison)
mcols(gr)$score <- c(1 - df$majiq$probability_changing,  df$leafcutter$p.adjust, df$junctionseq$padjust)

if (!all(as.logical(lapply(files, file.exists)))){
  stop("Input file not found.", call.=FALSE)
} else if (!file.exists(opt$annotation)) {
  stop("Annotation not found.", call.=FALSE)
}

message('Loading input')
files <- strsplit(opt$input, ',')[[1]]
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

gr_metadata <- mcols(gr)[queryHits(hits), c('comparison', 'method', 'score')]
gr_metadata$idx <- subjectHits(hits)
gr_metadata_agg <- aggregate(. ~ idx, gr_metadata, paste, collapse = ';')

message('Matching SJ to introns')

introns_with_match <- merge(
  as.data.frame(introns),
  gr_metadata_agg,
  by.x = "row.names",
  by.y = "idx",
  no.dups = T
)
introns_with_match <- subset(introns_with_match, select = -Row.names)
write_csv(introns_with_match, opt$output)
#introns_wo_match <- data[setdiff(seq_along(data), unique(queryHits(hits))), ]
#write.csv(introns_wo_match, sub('csv', opt$output, 'nomatch.csv'))
