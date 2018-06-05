# Title     : TODO
# Objective : TODO
# Created by: tbrittoborges
# Created on: 05/06/2018

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

save.cell.plot <- function(x, title){
  x <- subset(x, pvalCutOff <= 0.05)
  x <- head(x[order(-x$LogEnrich),], 30)

  pdf(paste0(title, "cell_plot.pdf"), 7, 5)
  cell.plot(x = setNames(x$LogEnrich, x$Term),
            cells = x$log2FoldChange,
            main = title)
  dev.off()

  pdf(paste0(title, "sym_plot.pdf"), 7, 5)
  sym.plot(x = setNames(x$LogEnrich, x$Term),
           cells = x$log2FoldChange,
           x.annotated = x$Annotated,
           main = title)
  dev.off()
  }


cell.plot.proc <- function(DESeq.res, ontology = 'BP') {
  DESeq.res <- as.data.frame(DESeq.res)
  DESeq.res <- na.omit(DESeq.res)
  is.de <- DESeq.res$padj<0.05

  all.genes <- DESeq.res$GeneID
  interesting_genes <- DESeq.res[DESeq.res$padj<0.05, 'GeneID']
  geneList <- factor(as.integer (all.genes %in% interesting_genes))
  names (geneList) <- all.genes

  g <- new(
    "topGOdata",
    ontology = ontology,
    allGenes =  geneList,
    mapping =  "org.Mm.eg.db",
    geneSelectionFun = function (allScore) {
      allScore <= 0.05
    },
    annotationFun = annFUN.org,
    ID = "Ensembl"
  )
  t <- new("elimCount",
           testStatistic = GOFisherTest,
           name = "Fisher test")
  s <- getSigGroups(g, t)
  r <- GenTable(g,
                pvalCutOff = s,
                topNodes = length(g@graph@nodes))
  r$pvalCutOff <- as.numeric(
    str_replace_all(r$pvalCutOff, "[^0-9e\\-\\.]*", ""))
  r$LogEnrich <- log2(r$Significant / r$Expected)
  ga <- genesInTerm(g)
  ga <- ga[r$GO.ID]
  names(ga) <- NULL
  r$GenesAnnotated <- ga
  xs <- DESeq.res[, c("padj", "log2FoldChange")]
  xs <- subset(xs, padj < 0.05)
  r$GenesSignificant <-
    lapply(r$GenesAnnotated, intersect, rownames(xs)) # extract genes
  ei.rows <- mclapply(r$GenesSignificant, function (y) {
    if (length(y))
      as.list(xs[y, , drop = FALSE])
    else
      as.list(rep(NA_real_, length(xs)))
  }, mc.cores = 10)
  ei <- mclapply(names(xs), function(z) {
    lapply(ei.rows, "[[", z)
  }, mc.cores = 10)
  ei <-
    structure(
      ei,
      names = names(xs),
      row.names = seq(nrow(r)),
      class = "data.frame"
    )
  row.names(ei) <- NULL
  r <-
    data.frame(r, ei, stringsAsFactors = FALSE, check.names = FALSE)
  return(r)
}

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
  res$GeneID <- row.names(res)

  # dir.create(base, recursive = T, showWarnings = F)

  # is.de <- as.numeric(
  #   x$padj < 0.05 & (x$log2FoldChange > 0.05 | x$log2FoldChange < 0.05))
  is.de <- as.numeric(res$padj<0.05)

  # glMDPlot(
  #   res,
  #   counts = dge$counts,
  #   status = is.de,
  #   main = title,
  #   groups = group,
  #   launch = F,
  #   anno = genes,
  #   folder = base)

  write.table(res,
              paste0(title,  "_de_genes.csv"),
              sep = ",")

  # bp.go.terms <- cell.plot.proc(res)
  # mf.go.terms <- cell.plot.proc(res, ontology = 'MF')

  # save.cell.plot(bp.go.terms,
  #                title= paste('BP_', title))
  # save.cell.plot(mf.go.terms,
  #                title = 'MF_' + title)
}

contrasts <- read.table(
  '/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/workflow_for_Thiago/contrasts.txt',
  header = T
)


apply(contrasts, 1, analysis.for)

# DEseq plots
#
# alpha <- 0.05
# cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
# plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
#      main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
#      pch=20, cex=0.6)
# abline(v=0)
# abline(v=c(-1,1), col="brown")
# abline(h=-log10(alpha), col="brown")
#
# gn.selected <- abs(res$log2FoldChange) > 2 & res$padj < alpha
# text(res$log2FoldChange[gn.selected],
#      -log10(res$padj)[gn.selected],
#      lab=rownames(res)[gn.selected ], cex=0.4)

# par(mfrow=c(2,1),cex.lab=0.7)
# epsilon=1
# col.pheno.selected <- group
# dds.norm <-  estimateSizeFactors(dds)
# boxplot(log2(counts(dds.norm)+epsilon),  col=col.pheno.selected,
#         las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
# boxplot(log2(counts(dds.norm, normalized=TRUE)+epsilon),  col=col.pheno.selected,
#         las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts")

# prop.null <- apply(dge, 2, function(x) 100*mean(x==0))
# barplot(prop.null, main="Percentage of null counts per sample",
#         horiz=TRUE, cex.names=0.5, las=1
#         , ylab='Samples', xlab='% of null counts')

# stats.per.sample <- data.frame(do.call(cbind, lapply(dge$counts, summary)))
# stats.per.sample


# epsilon <- 1 # pseudo-count to avoid problems with log(0)
# hist(as.matrix(log2(df + epsilon)), breaks=100, col="blue", border="white",
#      main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes",
#      las=1, cex.axis=0.7)

# boxplot(log2(df + epsilon),  pch=".",
#         horizontal=TRUE, cex.axis=0.5,
#         las=1, ylab="Samples", xlab="log2(Counts +1)")