#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result",,
    metavar = "character"
  ),
  make_option(
    c("-a", "--annotation"),
    type = "character",
    help = "Path to annotation (GTF/GFF format",
    metavar = "character"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

gr  <- GRanges(opt$input)
gtf <- rtracklayer::import.gff2(opt$annotation)

tx <- subset(gtf, type == 'transcript')
tx <- subsetByOverlaps(tx, gr)
gtf <- gtf[gtf$transcript_id %in% unique(tx$transcript_id)]

ex <- subset(gtf, type == 'exon')
names(ex) <- ex$transcript_id
ex_tx <- split(ex, names(ex))
introns <- psetdiff(unlist(range(ex_tx), use.names = FALSE), ex_tx)
tx_id <- names(introns)
introns <- unlist(introns)

txid_to_gene <- setNames(tx$gene_name, nm = tx$transcript_id)
txid_to_tx_name <- setNames(tx$cmp_ref, nm = tx$transcript_id)
txid_to_classcode <- setNames(tx$class_code, nm = tx$transcript_id)

mcols(introns)$gene_name <- txid_to_gene[introns$tx_id]
mcols(introns)$transcript_name <- txid_to_tx_name[introns$tx_id]
mcols(introns)$class_code <- txid_to_classcode[introns$tx_id]
names(introns) <- NULL

write.csv(as.data.frame(introns), paste0(opt$input, 'parsed.csv'))