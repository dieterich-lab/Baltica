#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(optparse)
})

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Differently spliced intron file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="Annotation in the GTF/GFF format", metavar="character")
  # make_option(c("-t", "--type"), type="character", default="acceptor",
  #             help="Extract donnor or acceptor [default \"%default\"]",
  #             metavar="character")

);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (!file.exists(opt$input)){
  print_help(opt_parser)
  stop("Pass a valid intron file to -i or --input", call.=FALSE)
}
if (!file.exists(opt$annotation)) {
  print_help(opt_parser)
  stop("Pass a valid annotation to -a or --annotation ", call.=FALSE)
}
# if (!opt$type %in% c('donor', 'acceptor')){
#   print_help(opt_parser)
#   stop("Valid choices for -t --type parameters are [donor, acceptor]", call.=FALSE)
# }

select.exons <- function(pos, exons, flank.dist=150) {
  dist <- distanceToNearest(pos, exons)

  if (any(is.na(subjectHits(acceptor.dist)))) {
    warning('While trying to find nearest feature, missing values were created')
  }

  e <- exons[subjectHits(dist), ]
  annotated <- mcols(dist)$distance == 0
  # for novel exons, use junction start + flank.dist nt for as the novel exon
  novel <- pos[queryHits(dist), ][!annotated, ]
  novel <- flank(novel, dist, start = FALSE)
  genes <- unique(mcols(e)$gene_id)

  mcols(novel) <- mcols(exons[!annotated])
  e[!annotated] <- novel

  names(e) <- mcols(e)$gene_id
  score(e) <- 0

  exons.in.genes <- exons[mcols(exons)$gene_id %in% genes]
  e.negative <- GenomicRanges::setdiff(exons.in.genes, acceptor.exons)

  list(
    positive=e,
    negative=e.negative
  )
}


intron <- GRanges(opt$input)
gtf <- rtracklayer::import.gff(opt$annotation)

exons <- gtf[gtf$type == 'exon']
exons <- exons[!duplicated(exons)]
# sequence.length <- quantile(width(exons), 0.5)

fivepss <- ifelse(strand(intron) == '+', end(intron), start(intron))

acceptor <- GRanges(seqnames = seqnames(intron),
                    IRanges(fivepss, width=1),
                    strand = strand(intron))

acceptor.exons <- select.exons(acceptor, exons)
rtracklayer::export(acceptor.exons['positive'], 'positive_exons.bed')
rtracklayer::export(acceptor.exons['negative'], 'negative_exons.bed')

threepss <- ifelse(strand(intron) == '+', start(intron), end(intron))
donor <- GRanges(seqnames = seqnames(intron),
                 IRanges(threepss, width=1),
                 strand = strand(intron))

donor.exons <- select.exons(donor, exons)
rtracklayer::export(donor.exons['positive'], 'positive_exons.bed')
rtracklayer::export(donor.exons['negative'], 'negative_exons.bed')





# old #####

# tst.dist <- distance(intron, exons[tst, ])
# tst.dist
#
# intron.foward <- intron[strand(intron) == "+"]
# intron.reverse <- intron[strand(intron) == "-"]
#
# donor.exon.foward <- findOverlaps(
#   resize(intron.foward, 1, fix='start'),
#   resize(exons, 1, fix='end'),
#   select = 'first')
#
# sum(is.na(donor.exon.foward))
#
# head(donor.exon.foward)
#
# tst <- split(
#   subjectHits(donor.exon.foward),
#   queryHits(donor.exon.foward))
#
# acceptor.exon.foward <- findOverlaps(
#   resize(intron.foward, 1, fix='end'),
#   resize(exons, 1, fix='start'))
#
# donor.exon.reverse <- findOverlaps(
#   resize(intron.reverse, 1, fix='end'),
#   resize(exons, 1, fix='start'))
#
# acceptor.exon.reverse <- findOverlaps(
#   resize(intron.reverse, 1, fix='start'),
#   resize(exons, 1, fix='end'))
#
# length(donor.exon.foward)
# length(donor.exon.reverse)
# length(acceptor.exon.foward)
# length(acceptor.exon.reverse)
#
# # donor
# donor <- exons[c(
#   subjectHits(donor.exon.foward),
#   subjectHits(donor.exon.reverse))]
#
# names(donor) <- c(
#   intron.foward[ queryHits(donor.exon.foward), 'junction'],
#   intron.reverse[ queryHits(donor.exon.reverse), 'junction'])
# score(donor) <- 0
#
# rtracklayer::export(donor, 'Baltica/leafcutter/donor_exon.bed')
#
# # acceptor
# acceptor <- exons[c(
#   subjectHits(acceptor.exon.foward),
#   subjectHits(acceptor.exon.reverse))]
#
# names(acceptor) <- c(
#   intron.foward[ queryHits(acceptor.exon.foward), 'junction'],
#   intron.reverse[ queryHits(acceptor.exon.reverse), 'junction'])
# score(acceptor) <- 0
#
# rtracklayer::export(acceptor, 'Baltica/leafcutter/acceptor_exon.bed')
#
#
# names(intron) <- intron$junction
# score(intron) <- 0
# rtracklayer::export(intron, 'Baltica/leafcutter/intron.bed')
#
#
