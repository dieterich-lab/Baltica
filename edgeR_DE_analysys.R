# Title     : edgeR differential expression analysis pipeline
#             with GO term analysys (gProfileR) and HTML
#             reports (Glimma)
# Objective : TODO
# Created by: tbrittoborges
# Created on: 24/05/2018

#!/usr/bin/env Rscrip
suppressPackageStartupMessages({
  library("stringr")
  library('ggplot2')
  library('knitr')
  library("edgeR")
  library('dplyr')
  library('RColorBrewer')
  library('pheatmap')
  library('DESeq2')
  library("regionReport")
  library('gProfileR')
  library("biomaRt")
  library("Glimma")
})

setwd('~/Desktop/TAC_DE/')
dir.create('DEanalysis/report/',
           showWarnings = FALSE,
           recursive = TRUE)

df <- read.csv(
  "gene_counts.csv",
  header = T,
  row.names = 1)

if (all(startsWith(rownames(df), 'ENSMUS'))){
  org <- 'mmusculus'
  library(Mus.musculus)
  org.lib <- Mus.musculus
} else if (all(startsWith(rownames(df), 'ENSG'))){
  org <- 'hsapiens'
  library(Homo.sapiens)
  org.lib <- Homo.sapiens
}
# else break



# sort the columns
df <- df[,c("Sham3h_1", "Sham3h_2", "Sham3h_3", "Tac3h_1", "Tac3h_2", "Tac3h_3",
      "Sham2d_1", "Sham2d_2", "Tac2d_1", "Tac2d_2",  "Tac2d_3", "Tac2d_4",
      "Tac2w_1",  "Tac2w_2",  "Tac2w_3")]
group <- factor(str_split(colnames(df), "_", simplify = T)[, 1])
dge <- DGEList(counts = df, group = group)

# remove genes with low CPM
dim(dge)
cpm <- cpm(dge)
keep.exprs <- rowSums(cpm>1)>=3
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)
dge <- calcNormFactors(dge, method="TMM")
# 10000000 / min(dge$samples$lib.size)

glMDSPlot(
  dge,
  labels=colnames(df),
  groups=group,
  launch=T,
  folder='DEanalysis/report/glMDSPlot')

design <- model.matrix(~0 + group)

dge <- estimateDisp(dge, design, robust = T)
plotBCV(dge)


cond.a <- colnames(design)[5]
cond.b <- colnames(design)[2]
cat('contrats is', cond.a, ' vs ', cond.b)
contr.matrix <- makeContrasts(
  paste(cond.a, cond.b, sep = '-'),
  levels = colnames(design)
)

res <- glmQLFTest(fit, contrast=contr.matrix)

is.de <- decideTestsDGE(res)
summary(is.de)

tr <- glmTreat(fit, contrast=contr.matrix, lfc=log2(2))
topTags(tr)

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

genes <- select(
  org.lib,
  keys = rownames(tr),
  columns = c("SYMBOL", "UNIPROT"),
  keytype = "ENSEMBL"
)
genes <- genes[!duplicated(genes$ENSEMBL),]
genes$UNIPROT <- paste0("<a href='https://www.uniprot.org/uniprot/",
                        genes$UNIPROT,
                        ".html'>",
                        genes$UNIPROT,
                        "</a>")

genes$ENSEMBL <- paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=",
                        genes$ENSEMBL,
                        "'>",
                        genes$ENSEMBL,
                        "</a>")


glMDPlot(tr,
         counts = dge$counts,
         coef = 1,
         status = is.de,
         main = colnames(tr)[1],
         groups = group,
         launch = T,
         anno = genes)

#
# subset(dt, dt == 1)
# genes <- select(org.lib,
#                 keys=rownames(dge),
#                 columns=c("SYMBOL", "DEFINITION"),
#                 keytype="ENSEMBL")
#
#
#
# #edge
# edgeReport(
#   dge,
#   lrt,
#   project = paste(unique(group), sep = " vs "),
#   intgroup = 'group',
#   browse = TRUE,
#   outdir = 'DEanalysis/edgeReport/')
#
# sig <- topTags(
#   lrt,
#   n = Inf,
#   p.value = 0.05,
#   sort.by = "logFC",
#   adjust.method = 'fdr')
# write.table(
#   sig,
#   "DEanalysis/edgeReport/sig_genes.txt",
#   sep = "\t")
#
# # plotSmear(lrt)
#
#
# sel <- row.names(sig)
# universe <- rowSums(cpm(dge) >= 1) >= 2
# universe <- as.vector(names(universe[universe]))
#
# goterm <- gprofiler(query = sel, organism = org, sort_by_structure = T, ordered_query = T, significant = T, exclude_iea = T, underrep = F, evcodes = F, region_query = F, max_p_value = 1, min_set_size = 3, max_set_size = 500, min_isect_size = 0, correction_method = "fdr", hier_filtering = "none", domain_size = "annotated", custom_bg = universe, numeric_ns = "", include_graph = F, src_filter = NULL)
# write.table(goterm, "DEanalysis/goterm.txt", sep = "\t")