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
# 44 negative 
pos = parse_suppa_sup(url, 4, event_id = TRUE, startRow = 3)
neg = parse_suppa_sup(url, 10, event_id = FALSE, 3)

neg <- neg %>%
    rownames_to_column('event_n') %>%
    dplyr::rename(chr=Chr, Gene=Gene_symbol)

reference <- rtracklayer::import.gff("/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.SIRV.gtf")

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
ch <- rtracklayer::import.chain('~/GRCh37_to_GRCh38_2.chain')

ref_ex_tx <- filter_multi_exon(reference)
ref_introns <- get_introns(ref_ex_tx)
ref_introns <- ref_introns[width(ref_introns) > 2, ]

neg_gr <- GRanges(
    seqnames = neg$chr, 
    ranges = IRanges(neg$Exon_start, neg$Exon_end))
seqlevelsStyle(neg_gr) <- 'ensembl'
neg_gr <- rtracklayer::liftOver(neg_gr, ch) 
neg_gr <- unlist(neg_gr)

neg_introns <- subsetByOverlaps(ref_introns, neg_gr)
# perc_olaps <- width(pintersect(ref_introns[from(neg_olaps)], neg_gr[to(neg_olaps)])) / width(neg_gr[to(neg_olaps)])
# neg_introns <- neg_introns[ perc_olaps == 1 ]
score(neg_introns) <- 0

pos <- pos %>% 
    rownames_to_column('event_n') %>% 
    mutate(
        # exon_coord = str_glue_data(., "{exon_start}-{exon_end}"),
        exon_start = NULL,
        exon_end = NULL,
        Event_id = NULL, 
        coord=str_glue('{coord1}|{coord2}'),
        coord1 = NULL,
        coord2 = NULL,
        exon_coord = NULL) %>% 
    separate_rows(coord, sep =  '\\|') %>% 
    mutate(coord = as.character(coord)) %>% 
    filter(coord != '')

gr <- GRanges(seqnames = pos$chr, ranges = IRanges(pos$coord), strand = pos$strand)
seqlevelsStyle(gr) <- 'ensembl'
# https://gist.github.com/tbrittoborges/f0bc89cd7192d0955679d5105fd39475
gr <- rtracklayer::liftOver(gr, ch) 
gr <- unlist(gr)

pairs <- findOverlapPairs(ref_introns, gr)
same_end <- end(pairs@first) - end(pairs@second) == -1
same_start <- start(pairs@first) - start(pairs@second) == 1
length(unique(pairs[same_start & same_end]@second))
pairs <- pairs[same_start & same_end]
gr <- unique(pairs@first)
mcols(gr)$score <- 1
names(gr) <- seq_along(gr)

gr <- c(neg_introns, gr)
export.bed(gr, '/prj/Niels_Gehring/baltica_benchmark/best_et_al/baltica/ortho.bed')
