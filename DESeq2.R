#!/usr/bin/env Rscrip
suppressPackageStartupMessages({
  library(stringr)
  library(DESeq2)
  library(IHW)
  library(ggplot2)
  library(clusterProfiler)
})

cts <- read.csv(
  "/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/gene_count_matrix.csv",
  header = T, 
  row.names = 1)
setwd('/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/DESeq2/')

vals <- as.vector(t(str_match(colnames(cts), '(.*)_(.*)') [, -1]))
condition <- as.factor(vals[seq(1, length(vals), 2)])
coldata <- data.frame(row.names=colnames(cts), 
                      condition=condition)
all(rownames(coldata) == colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ condition)

analysis <- function(contrast) {
  dds$condition <- relevel(dds$condition, ref=as.character(contrast[1]))
  dds.this <- DESeq(dds)
  
  coef <- grep(as.character(contrast[2]), resultsNames(dds.this))
  res.ihw <- results(
    dds.this, 
    filterFun=ihw, 
    contrast = c("condition", 
                 as.character(contrast[2]), 
                 as.character(contrast[1])), 
    lfcThreshold=.5, 
    altHypothesis="greaterAbs",
    alpha=0.05)
  res <- lfcShrink(
    dds.this, 
    res=res.ihw, 
    type='apeglm', 
    coef=coef)
  
  write.csv(as.data.frame(res[
    with(subset(res, padj < 0.05),
         order(abs(res$log2FoldChange))), ]),
    file=paste(
      as.character(contrast[2]), 
      as.character(contrast[1]),
      'DESeq2.csv', 
      collapse = '_'))
  # write.csv(
  #   summary(res),
  #   file=paste(
  #     as.character(contrast[2]), 
  #     as.character(contrast[1]),
  #     'DESeq2_summary.csv', 
  #     collapse = '_'))
    
  
  # from https://www.r-bloggers.com/deseq2-course-work/
  
  res.df <- as.data.frame(res)
  sig=ifelse(res.df$padj<0.05, "FDR<0.05", "Not Sig")
  p = ggplot2::ggplot(res.df, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(col = sig)) +
    ggplot2::scale_color_manual(values = c("red", "black")) +
    ggplot2::ggtitle("Volcano Plot of DESeq2 analysis" )
  top <- head(
    subset(res.df[order(abs(res.df$log2FoldChange), decreasing = TRUE), ], padj < 0.05),
    n=10)
  symbols <- bitr(
    row.names(top), 'ENSEMBL', 'SYMBOL', "org.Mm.eg.db", drop = FALSE)
  symbols$SYMBOL[is.na(symbols$SYMBOL)] <- symbols$ENSEMBL[is.na(symbols$SYMBOL)] 
  p + ggrepel::geom_text_repel(data=top, ggplot2::aes(label=symbols$SYMBOL))
  ggplot2::ggsave(paste0(
        as.character(contrast[2]), 
        as.character(contrast[1]),
        'DESeq2.png', collapse = '_'))
  
}
contrasts <- read.table('/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/workflow_for_Thiago/contrasts.txt', header = T)
apply(contrasts, 1, analysis)