# Title     : TODO
# Objective : TODO
# Created by: tbrittoborges
# Created on: 05/06/2018

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
design <- model.matrix(~ 0 + group)
v <- voom(x, design, plot = F)
vfit <- lmFit(v, design)

analysis <- function(variables) {
  cond.name <- paste0(
    'Limma_', variables[1], '_vs_', variables[2])
  a <- paste0('group', variables[1])
  b <- paste0('group', variables[2])
  c <- paste(a, b, sep = '-')
  cat(c,  "\n")

  contr.matrix <- makeContrasts(
    contrasts = c, levels = factor(colnames(design)))

  vfit2 <- contrasts.fit(
    vfit, contrasts = contr.matrix)
  efit <- eBayes(vfit2, proportion = 0.05)
  res <- topTable(efit, coef=1, adjust="BH", number = Inf)
  res$source <- cond.name
  write.table(res,
              paste(cond.name,  "_de_genes.csv", sep = ""),
              sep = ",")
}

contrasts <- read.table('/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/workflow_for_Thiago/contrasts.txt', header = T)
apply(contrasts, 1, analysis)