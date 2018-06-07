suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(topGO)
  library(CellPlot)
  library(stringr)
  library(Glimma)
})
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
setwd('/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/')

df <- read.csv("gene_count_matrix.csv",
               header = T,
               row.names = 1)

if (all(startsWith(rownames(df), 'ENSMUS'))) {
  org <- 'mmusculus'
  library(Mus.musculus)
  org.lib <- Mus.musculus
} else if (all(startsWith(rownames(df), 'ENSG'))) {
  org <- 'hsapiens'
  library(Homo.sapiens)
  org.lib <- Homo.sapiens
}

vals <- as.vector(t(str_match(colnames(df), '(.*)_(.*)') [, -1]))
group <- as.factor(vals[seq(1, length(vals), 2)])
dge <- DGEList(counts = df, group = group)
keep.exprs <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[keep.exprs, , keep.lib.sizes = FALSE]

genes <- select(
  org.lib,
  keys = rownames(dge),
  columns = c("SYMBOL", "UNIPROT"),
  keytype = "ENSEMBL"
)
genes <- genes[!duplicated(genes$ENSEMBL), ]
row.names(genes) <- genes$ENSEMBL
genes$UNIPROT <- paste0(
  "<a href='https://www.uniprot.org/uniprot/",
  genes$UNIPROT,
  ".html'>",
  genes$UNIPROT,
  "</a>"
)

genes$ENSEMBL <-
  paste0(
    "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=",
    genes$ENSEMBL,
    "'>",
    genes$ENSEMBL,
    "</a>"
  )
genes$GeneID <- rownames(genes)

dds <- DESeqDataSetFromMatrix(
  as.matrix(dge), 
  DataFrame(group), ~ group)
dds <- DESeq(dds)

analysis.for <- function(x) {
  title <- paste0('DESeq2_', x[1], '_vs_', x[2])
  cat(title, "\n")
  res <- results(
    dds,
    contrast = c('group', x[1], x[2]),
    alpha = 0.05,
    pAdjustMethod = "BH"
  )
  res$source <- gsub('DESeq2_', '', title)
  base <- paste0('DESeq2_', title, '/', title)
  write.table(
    as.data.frame(res),
    paste(title,  "_de_genes.csv", sep = ""),
    sep = ",")

  # https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#exporting-results-to-glimma
  res.df <- as.data.frame(res)
  res.df$log10MeanNormCount <- log10(res.df$baseMean)
  idx <- rowSums(counts(dds)) > 0
  res.df <- res.df[idx,]
  
  res.df$padj[is.na(res.df$padj)] <- 1
  is.de <- as.numeric(res.df$padj<0.05)
  res.df$GeneID <- rownames(res.df)
    
  glMDPlot(
    res.df,
    xval="log10MeanNormCount",
    yval="log2FoldChange",
    counts=counts(dds)[idx,],
    anno=genes,
    groups=group,
    samples=colnames(dds),
    status=is.de,
    folder = title)
}
analysis.for(c("Sham_3h_tac", "TAC_3h"))
contrasts <- read.table(
  '/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/workflow_for_Thiago/contrasts.txt',
  header = T
)
apply(contrasts, 1, analysis.for)


# res <- results(
#   dds,
#   contrast = c('group', "Sham_3h_tac",  "TAC_3h"),
#   alpha = 0.05,
#   pAdjustMethod = "BH"
# )
# res.df <- as.data.frame(res)
# res.df$log10MeanNormCount <- log10(res.df$baseMean)
# idx <- rowSums(counts(dds)) > 0
# res.df <- res.df[idx,]
# res.df$padj[is.na(res.df$padj)] <- 1
# glMDPlot(
  # res.df,
  # xval="log10MeanNormCount",
  # yval="log2FoldChange",
  # counts=counts(dds),
  # anno=genes[rownames(dds), ],
  # groups=group,
  # samples=colnames(dds),
  # status=res.df$padj < 0.1,
  # folder = title)