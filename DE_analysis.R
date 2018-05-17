#!/usr/bin/env Rscrip
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('knitr'))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('RColorBrewer'))
suppressPackageStartupMessages(library('pheatmap'))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library("regionReport"))
suppressPackageStartupMessages(library('gProfileR'))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("Glimma"))

dir.create('DEanalysis/limmaReport/',
           showWarnings = FALSE,
           recursive = TRUE)

plot.after.normal <- function() {
  # verify genecount before and after
  nsamples <- ncol(dge)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow = c(1, 2))
  plot(
    density(lcpm[, 1]),
    col = col[1],
    lwd = 2,
    ylim = c(0, 0.21),
    las = 2,
    main = "",
    xlab = ""
  )
  title(main = "A. Raw data", xlab = "Log-cpm")
  abline(v = 0, lty = 3)
  for (i in 2:nsamples) {
    den <- density(lcpm[, i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright",
         legend = group,
         text.col = col,
         bty = "n")
  lcpm <- cpm(dge, log = TRUE)
  plot(
    density(lcpm[, 1]),
    col = col[1],
    lwd = 2,
    ylim = c(0, 0.21),
    las = 2,
    main = "",
    xlab = ""
  )
  title(main = "B. Filtered data", xlab = "Log-cpm")
  abline(v = 0, lty = 3)
  for (i in 2:nsamples) {
    den <- density(lcpm[, i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright",
         legend = group,
         text.col = col,
         bty = "n")

}

df <- read.csv("DEanalysis/gene_counts.csv",
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

group <- factor(str_split(colnames(df), "_", simplify = T)[, 1])
x <- DGEList(counts = df, group = group)
par(mfrow = c(1, 1))
# plotMDS(dge)

# remove genes with low CPM
cpm <- cpm(x)
lcpm <- cpm(x, log = TRUE)
glMDSPlot(
  lcpm,
  labels = colnames(df),
  groups = group,
  launch = F,
  folder = 'DEanalysis/limmaReport/glMDSPlot'
)
# plot.after.normal()
dim(x)
keep.exprs <- rowSums(cpm > 1) >= 3
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
# plot.after.normal()
dim(x)

x <- calcNormFactors(x, method = 'TMM')
design <- model.matrix(~ 0 + group)

cond.a <- colnames(design)[1]
cond.b <- colnames(design)[2]
contr.matrix <- makeContrasts(
  paste(cond.a, cond.b, sep = '-'),
  levels = colnames(design))

v <- voom(x, design, plot = F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(
  vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
# plotSA(efit)

summary(decideTests(efit))
tfit <- treat(vfit, lfc = .5)
dt <- decideTests(tfit)
summary(dt)

genes <- select(
  org.lib,
  keys = rownames(x),
  columns = c("SYMBOL", "UNIPROT"),
  keytype = "ENSEMBL"
)
genes <- genes[!duplicated(genes$ENSEMBL),]

genes$UNIPROT <- paste0("<a href='https://www.uniprot.org/uniprot/",
       genes$UNIPROT,
       ".html'>",
       genes$UNIPROT,
       "</a>")

genes$ENSEMBL <- paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", genes$ENSEMBL,

                        "'>",
                        genes$ENSEMBL,
                        "</a>")

glMDPlot(
  tfit,
  anno = genes,
  coef = 1,
  status = dt,
  main = colnames(tfit)[1],
  counts = x$counts,
  groups = group,
  launch = T,
  folder = 'DEanalysis/limmaReport/glMDPlot'
)

sel <- row.names(subset(dt, dt != 1))
universe <- as.vector(rownames(x))

goterm <- gprofiler(
  query = sel,
  custom_bg = universe,
  organism = org,
  sort_by_structure = T,
  ordered_query = T,
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
  numeric_ns = "",
  include_graph = T,
  src_filter = NULL
)
write.table(goterm,
            "DEanalysis/goterm.txt", sep = "\t")

pheatmap(x[sel, ], color=rev(brewer.pal(6, "RdBu")), scale="row",
         kmeans_k=NA, show_rownames=F, show_colnames=T,
         main="Sham vs Tac (3h))", cluster_rows=TRUE,
         cluster_cols=FALSE, clustering_distance_rows="correlation")


