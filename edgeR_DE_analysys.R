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
  library('gProfileR')
  library("biomaRt")
  library("Glimma")
})

setwd('/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/')

df <- read.csv(
  "gene_count_matrix.csv",
  header = T,
  row.names = 1)

vals <- as.vector(t(str_match(colnames(df), '(.*)_(.*)') [, -1]))
group <- as.factor(vals[seq(1, length(vals), 2)])
x <- DGEList(counts = df, group = group)
# remove genes with low CPM
cpm <- cpm(x)
keep.exprs <- rowSums(cpm > 1) >= 3
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
x <- calcNormFactors(x, method = 'TMM')
x <- estimateDisp(x, design, robust=TRUE)

batch <- '1 2 2 1 1 1 1 1 1 1 3 2 2 2 1 1 2 2 2 1 1 1 2 2 2 1 1 1 1 1'
batch <- factor(str_split(batch, ' ', simplify = T))



glMDSPlot(
  dge,
  labels=colnames(df),
  groups=group,
  launch=T,
  folder='glMDSPlot')

genes <- select(
  org.lib,
  keys = rownames(dge),
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
design <- model.matrix(~group)
dge <- estimateDisp(dge, design, robust = T)
dge <- calcNormFactors(dge)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}


analysis <- function(variables) {
    title <- paste0('Edger_', variables[1], '_vs_', variables[2])
    cat(title, "\n")
    dir.create(title,
               showWarnings = FALSE,
               recursive = TRUE)
  res <- exactTest(x, variables)
  is.de <- decideTestsDGE(res, p.value=.05)
  glMDPlot(res,
           counts = dge$counts,
           status = is.de,
           main = title,
           groups = group,
           launch = F,
           anno = genes,
           folder=title)
  write.table(is.de,
              paste(title,  "_de_genes.csv", sep = ""),
              sep = ",")

  top <- topTags(res, sort.by="logFC", n=30)
  pdf(paste(base, "_pheatmap.pdf", sep = ''), 7, 5)
  data_norm <- t(
    apply(
      dge$counts[rownames(top) ,], 1, cal_z_score))

  pheatmap(
    data_norm,
    color=rev(brewer.pal(6, "RdBu")),
    scale="row",
    kmeans_k=NA,
    show_rownames=T,
    show_colnames=T,
    cluster_rows=TRUE,
    cluster_cols=T,
    main = title,
    clustering_distance_rows="correlation")
  dev.off()

  up <- row.names(subset(is.de, is.de == 1))
  down <- row.names(subset(is.de, is.de == -1))
  universe <- as.vector(row.names(dge))

  goterm.up <- gprofiler(
    query = up,
    organism = org,
    sort_by_structure = T,
    ordered_query = F,  # TODO
    significant = T,
    exclude_iea = T,
    underrep = F,
    evcodes = F,
    region_query = F,
    max_p_value = 1,
    min_set_size = 3,
    max_set_size = 500,
    min_isect_size = 0,
    correction_method = "fdr",
    hier_filtering = "none",
    domain_size = "annotated",
    custom_bg = universe,
    numeric_ns = "",
    include_graph = F,
    src_filter = NULL)
  write.table(goterm.up,
              paste(base, "_upGoTerms.csv", sep = ""),
              sep = ',')

  goterm.down <- gprofiler(
    query = down,
    organism = org,
    sort_by_structure = T,
    ordered_query = F,  # TODO
    significant = T,
    exclude_iea = T,
    underrep = F,
    evcodes = F,
    region_query = F,
    max_p_value = 1,
    min_set_size = 3,
    max_set_size = 500,
    min_isect_size = 0,
    correction_method = "fdr",
    hier_filtering = "none",
    domain_size = "annotated",
    custom_bg = universe,
    numeric_ns = "",
    include_graph = F,
    src_filter = NULL)
  write.table(goterm.down,
              paste(base, "_downGoTerms.csv", sep = ""),
              sep = ',')
}
contrasts <- read.table('/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/workflow_for_Thiago/contrasts.txt', header = T)
apply(contrasts, 1, analysis)
