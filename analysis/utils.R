#' Compute and filter hits based on the fraction of overlap
#'
#' @param query first set of range
#' @param subject second set of ranges
#' @param cutoff for keeping the ranges
#' @return overlapping ranges given the contrain
#' @export
filter_hits_by_fraction <- function(query, subject, cutoff = 0.99) {
  stopifnot(is(query, "GRanges"))
  stopifnot(is(subject, "GRanges"))
  hits <- findOverlaps(query, subject)
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  overlap <- width(pintersect(query, subject)) / pmin(width(query), width(subject))
  hits[overlap >= cutoff]
}

#' Compute and filter hits based on the difference in the genomic start and end
#'
#' @param query first set of range
#' @param subject second set of ranges
#' @param max_start, max_end absolute max difference at the start and end coordinates, respectively
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
}


#' Compute a set of introns from GTF files, tested with the ones produced by StringTie
#'
#' @param gtf_path path to to the gtf file
#' @param read_gtf function to read the gtf_file
#' @return a GRange that obj with the introns named by their parent trancripts
#' @export
get_introns <- function(gtf_path, read_gtf = rtracklayer::import.gff2) {
  gtf <- import.gff2(gtf_path)
  tx <- subset(gtf, type == 'transcript')
  ex <- subset(gtf, type == 'exon')
  multi_ex <- table(ex$transcript_id) > 1
  ex <- subset(ex, mcols(ex)$transcript_id %in% names(multi_ex[multi_ex]))
  names(ex) <- ex$transcript_id
  ex_tx <- split(ex, names(ex))
  introns <- psetdiff(unlist(range(ex_tx), use.names = FALSE), ex_tx)
  introns_names <- names(introns)
  introns_len <- lapply(introns, length)
  mcols(introns)['transcript_id'] <- rep(introns_names, introns_len)
  introns <- unlist(introns)
  introns
}
