# Title     : Use TopGo gene ontology with DESeq2
# Created by: tbrittoborges
# Created on: 10/07/2018

library(CellPlot)
library(topGO)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(edgeR)
library(biomaRt)

setwd('/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/CellPlot_BP/')
files <- Sys.glob(
  '/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/DESeq2/*.csv')

# data for the background
gene_counts <- read.csv(
  '/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/gene_count_matrix.csv',
  header = 1,
  row.names = 1
)
cpm <- cpm(gene_counts)
dim(cpm)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

get.title <- function(file.name){
  title <- str_split(basename(file.name), ' ')[[1]]
  title <- paste( title[1:2], collapse = '_vs_')
  title
}
go.list <- list()
cell.plot.input <- function(f){
  # foreground
  x <- read.csv(f, row.names = 1)
  title <- get.title(f)
  fg.ds <- subset(x, x$padj < 0.05)
  foreground <- row.names(fg.ds)
  # background
  split.title <- strsplit(title, '_vs_')[[1]]
  background.colnames <- gsub('_', '', colnames(cpm))
  sel.cols <- lapply(split.title, grep, background.colnames)
  background <- names(rowSums(cpm[, unlist(sel.cols)]) > 3)

  # topgo pipeline from cellplot modified by Alexey
  inUniverse <- background %in% background
  inSelection <- background %in% foreground
  alg <- factor(as.integer(inSelection[inUniverse]))
  names(alg) <- background
  g <- new(
    "topGOdata",
    ontology = "BP",
    description = title,
    allGenes = alg,
    mapping = "org.Mm.eg.db",
    annotationFun = annFUN.org,
    ID = "Ensembl"
  )
  t <- new("elimCount",
           testStatistic = GOFisherTest,
           name = "Fisher test") # test definition
  s <- getSigGroups(g, t) # run F-test
  r <- GenTable(g,
                pvalCutOff = s,
                topNodes = length(g@graph@nodes)) # return data.frame

  r$pvalCutOff <- as.numeric(str_replace_all(r$pvalCutOff, "[^0-9e\\-\\.]*", ""))
  r$LogEnrich <- log2(r$Significant / r$Expected)
  write.csv(r, file=paste(title, 'topGO.csv', collapse = '_'))
  ga <- genesInTerm(g) # GenesAnnotated | list of genes per go-terms
  ga <- ga[r$GO.ID] # eliminate missing terms
  names(ga) <- NULL
  r$GenesAnnotated <- ga

  r$GenesSignificant <-
    lapply(r$GenesAnnotated, intersect, foreground) # extract genes

  r$padj <- lapply(r$GenesSignificant, function (y) {
    pval.match <- match(y, row.names(fg.ds))
    fg.ds$padj[pval.match]
  })

  r$lfc <- lapply(r$GenesSignificant, function (y) {
    pval.match <- match(y, row.names(fg.ds))
    fg.ds$log2FoldChange[pval.match]
  })
  r
}

sfbGO<-lapply(files, cell.plot.input)
names(sfbGO)<-lapply(
  files, get.title)

for (n in names(sfbGO)) {
  cat(n, '\n')
  x <- subset(sfbGO[[n]], pvalCutOff <= 0.05 & Significant > 30)
  x <- head(x[order(-x$LogEnrich),], n=30)

  try({
    png(filename = paste(n, 'cellplot.png', sep = '_'))
    cell.plot(x = setNames(x$LogEnrich, x$Term),
              cells = x$lfc,
              main = n,
              x.mar = c(.2, .1),
              key.n = 7,
              y.mar = c(.05, .1),
              cex = 1.6,
              cell.outer = 3,
              bar.scale = .5,
              space = .2)
    dev.off()
  })

  try({
    png(filename = paste(n, 'symplot.png', sep = '_'))
    sym.plot(x = setNames(x$LogEnrich, x$Term),
             cells = x$lfc,
             x.annotated = x$Annotated,
             main = n,
             y.mar = c(.05, .1),
             key.n = 7,
             # cex = 1.6,
             axis.cex = .8,
             group.cex = .7)
    dev.off()
  })

  try({
    x$up <- lapply(Map(setNames, x$lfc, x$GenesSignificant), function (i) { i[i>0] })
    x$dwn <- lapply(Map(setNames, x$lfc, x$GenesSignificant), function (i) { i[i<0] })
    png(filename = paste(n, 'arcplot.png',  sep = '_'))
    arc.plot(x = setNames(x$LogEnrich, x$Term),
             up.list = x$up,
             down.list = x$dwn,
             x.mar = c(.9, .5),
             main = n)
    dev.off()
  })
}

y <- lapply(sfbGO, function (xx) {
  xx$Upregulated <- unlist(lapply(xx$lfc, function (z) sum(z>0)))
  xx$Downregulated <- unlist(lapply(xx$lfc, function (z) sum(z<0)))
  xx
})

yterms <- unique(unlist(lapply(y, function(x){
  x <- subset(x, pvalCutOff <= 0.05)
  x <- x[order(x$LogEnrich),]
  head(x, 9)$GO.ID
})))

png(filename = 'gohistogram.png', width = 20, height = 20, res=300, units = 'cm')
go.histogram(y,
             go.alpha.term = "pvalCutOff",
             gene.alpha.term = "padj",
             logfc.term = "lfc",
             min.genes = 5,
             max.genes = 1e10,
             go.selection = yterms,
             show.ttest = T,
             main = "GO enrichment\n",
             axis.cex = 1,
             lab.cex = 1.5,
             main.cex = 1.5
             )
dev.off()

y2 <- lapply(sfbGO[c("TAC_3h_vs_Sham_3h_tac", "TAC_2d_vs_TAC_3h",  "TAC_2wks_vs_TAC_2d", "TAC_2wks_vs_Sham_2wks_tac" )], function (xx) {
  xx$Upregulated <- unlist(lapply(xx$lfc, function (z) sum(z>0)))
  xx$Downregulated <- unlist(lapply(xx$lfc, function (z) sum(z<0)))
  xx
})

yterms2 <- unique(unlist(lapply(y, function(x){
  x <- subset(x, pvalCutOff <= 0.05)
  x <- x[order(x$LogEnrich),]
  head(x, 9)$GO.ID
})))

png(filename = 'gohistogramTAC.png', width = 20, height = 20, res=300, units = 'cm')
go.histogram(y2,
             go.alpha.term = "pvalCutOff",
             gene.alpha.term = "padj",
             logfc.term = "lfc",
             min.genes = 5,
             max.genes = 1e10,
             go.selection = yterms,
             show.ttest = T,
             main = "GO enrichment\n",
             axis.cex = 1,
             lab.cex = 1.5,
             main.cex = 1.5
             # lab.cex = c('a', 'b', 'c', 'd','e','f')
)
dev.off()
