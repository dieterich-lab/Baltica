library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)
library(GenomicRanges)

parse_suppa_sup = function(url, sheet, event_id=TRUE, startRow = 3){
    x = openxlsx::read.xlsx(url, sheet = sheet, startRow = startRow)
    if (isTRUE(event_id)) {
        event_id = stringr::str_split(x$Event_id, ';|:', simplify = TRUE)
        colnames(event_id) = c('gene_id', 'type', 'chr', 'coord1', 'coord2', 'strand')
        event_id = as.data.frame(event_id)
        event_id$strand = NULL
        event_id$chr = NULL
        x = dplyr::bind_cols(x, event_id)
    }
    return(x)
}

url='https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-018-1417-1/MediaObjects/13059_2018_1417_MOESM2_ESM.xlsx'
# from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1417-1#MOESM2
# RT-PCR validated:
# 83 positive 
# 45 negative 
pos = parse_suppa_sup(url, 4, event_id = TRUE, startRow = 3)
neg = parse_suppa_sup(url, 10, event_id = FALSE, 3)

neg <- neg %>%
    rownames_to_column('event_n') %>%
    dplyr::rename(chr=Chr, Gene=Gene_symbol)

reference <- rtracklayer::import.gff("/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.SIRV.gtf")
reference <- subset(reference, source == 'ensembl_havana')

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

get_introns <- function(ex_tx) {
    stopifnot(is(ex_tx, "List"))

    introns <- psetdiff(unlist(range(ex_tx), use.names = FALSE), ex_tx)
    introns <- unlist(introns)
    introns
}

# see here how to process this 
ch <- rtracklayer::import.chain('~/GRCh37_to_GRCh38_2.chain')

ref_exons <- subset(reference, type == 'exon')
ref_exons <- unique(ref_exons)
genes <- subset(reference, type == 'gene')
genes <- as.data.frame(genes)[, c('strand', 'gene_name')]

neg <- left_join(neg, genes, by=c("Gene"="gene_name"))

neg_gr <- GRanges(
    seqnames = neg$chr, 
    ranges = IRanges(neg$Exon_start, neg$Exon_end), 
    strand = neg$strand)
score(neg_gr) <- 0
names(neg_gr) <- neg$Gene_symbol

pos_gr <- GRanges(
    seqnames = pos$chr,
    ranges = IRanges(pos$exon_start, pos$exon_end),
    strand = pos$strand)
score(pos_gr) <- 1
names(pos_gr) <- pos$Gene

gr <- c(pos_gr, neg_gr)
table(gr$score)

seqlevelsStyle(gr) <- 'ensembl'
gr2 <- gr
gr <- unlist(rtracklayer::liftOver(gr, ch))
strand(gr[strand(gr) == '*']) <- c('-', '+', '+', '-')

pre <- precede(gr, ref_exons)
fol <- follow(gr, ref_exons)
introns <- GRanges()
for(i in seq_along(gr)){
    both <- c(ref_exons[pre][i], ref_exons[fol][i])
    introns <- append(introns, setdiff(range(both), both))
}
score(introns) <- score(gr)
export.bed(introns, '/prj/Niels_Gehring/baltica_benchmark/best_et_al/baltica/ortho2.bed')
