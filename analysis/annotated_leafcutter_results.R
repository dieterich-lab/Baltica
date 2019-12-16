#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(GenomicRanges)
  library(openxlsx)
  library(yaml)
  library(optparse)
})

option_list = list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Leafcutter result parsed",
    metavar = "character",
  ),
  make_option(
    c("-a", "--annotation"),
    type = "character",
    default = NULL,
    help = "Annotation in the GTF/GFF format",
    metavar = "character",
  )
)

opt = parse_args(OptionParser(option_list=option_list))


if (!file.exists(opt$input)) {
  stop("Input file not found.", call.=FALSE)
} else if (!file.exists(opt$annotation)) {
  stop("Annotation not found.", call.=FALSE)
}

df <-  read_csv(
  opt$input,
  col_types = c(
    "chr" = "c",
    "genes" = "c",
    "contrast" = "c",
    "intron"  = "c",
    "strand" = "c",
    "cluster" = "c",
    "status"  = "c",
    .default = "d"
  )
)

# TODO pval and deltapsi as options

df <- df %>%
  filter(p.adjust < 0.05,  abs(deltapsi) > .1)


df[!df$strand %in% c('+', '-'), "strand"] <- '*'

gtf <- rtracklayer::import.gff(opt$annotation)
genes <- gtf[gtf$type == "gene"]
genes.list <- split(genes, genes$gene_id)

# find introns
# exons <- subset(gtf, type == "exon")
# exons.list <- GenomicRanges::split(exons, exons$gene_id)
# introns.list <- lapply(exons.list, gaps)
# introns <- unlist(as(introns.list, "GRangesList"))
# introns <- subsetByOverlaps(introns, genes, type='within')

# # ...
# hits <- findOverlaps(gr, introns)
# overlaps <- pintersect(gr[queryHits(hits)], introns[subjectHits(hits)])
# percentOverlap <- width(overlaps) / width(introns[subjectHits(hits)])
# # tst <- split(subjectHits(hits), queryHits(hits))
# gr['annotated'] <- percentOverlap == 1 # TODO ?
# # gr['canonical'] <- PSI ref column >  percentOverlap == 1 # TODO ?


gr <- GRanges(df[, c("chr", "start", "end", "strand")])
gr_on_gtf <- mergeByOverlaps(gr, genes)

df$junction <- as.character(str_glue_data(df, "{chr}:{start}-{end}:{strand}"))
gr_on_gtf$gr <- as.character(gr_on_gtf$gr)
gr_on_gtf <- as.data.frame(gr_on_gtf)
gr_on_gtf <-
  dplyr::select(gr_on_gtf, c("gr", "gene_id", "gene_name", "gene_biotype"))
df <- left_join(df, gr_on_gtf, by = c("junction" = "gr"))
df <- distinct(df)

up <- df %>%
  filter(p.adjust < 0.05, deltapsi > .1) %>%
  group_by(junction) %>%
  summarise(up = n_distinct(contrast))

do <- df %>%
  filter(p.adjust < 0.05, deltapsi < -.1) %>%
  group_by(junction) %>%
  summarise(down = n_distinct(contrast))

df <- left_join(df, up, by = "junction")
df <- left_join(df, do, by = "junction")
df$up <- tidyr::replace_na(df$up, 0)
df$down <- tidyr::replace_na(df$down, 0)
df["overlaps"] <- df$up - df$down

out.path <- dirname(opt$input)

message('Writting the annotated output to ', file.path(out.path, 'leafcutter_junctions_annotated.csv'))


write_csv(distinct(df), file.path(out.path, "leafcutter_junctions_annotated.csv"))

wb <-  createWorkbook()
for (cont in unique(df$contrast)) {
  addWorksheet(wb, cont)
  tmp <- df %>% filter(contrast == cont)
  tmp <- tmp[colSums(!is.na(tmp)) > 0]
  writeDataTable(wb, tmp, sheet = cont)
}

col_desc <- read.csv(
    text = "contrast,pair-wise comparison
    intron,intron identifier (refer to this contrast)
    logef,likelihood ratio: null model vs alt model
    deltapsi,PSIref - PSIalt
    chr,chromosome or contig
    start,genomic start position
    end,genomic end position
    strand,RNA strand
    cluster,unique identifier for the cluster (refer to this contrast)
    status,OK or reason not to tested
    loglr,log effect size
    df,degrees of freedom (number of introns within a cluster - 1)
    p,p-value
    p.adjust,adjusted p-value (method = FDR)
    genes,gene name
    junction,junction on the format chr:start-end:strand
    gene_id,gene identifier from Ensembl
    gene_name,gene symbol
    gene_biotype,gene biotype from Ensembl
    up,junction upregulated `up` times (all contrasts)
    down,junction downregulated `down` times (all contrasts)
    overlaps,up - down
    columns with condition name,PSI value (fitted junction usage)"
  )


addWorksheet(wb,  'col_desc')
writeDataTable(wb, sheet = 'col_desc', as.data.frame(col_desc))

session.info <- capture.output(sessionInfo())
addWorksheet(wb, 'session info')
writeDataTable(wb, sheet = 'session info', as.data.frame(session.info))

saveWorkbook(wb,  file.path(out.path, "DS_results_annotated.xlsx" ), overwrite=T)

