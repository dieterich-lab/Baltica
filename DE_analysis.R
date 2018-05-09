# ls /Volumes/prj/Niels_Gehring/newData_March_2018/workflow/mapping/*_STARmapping/ReadsPerGene.out.tab | grep -v 'mate'

# paste  ls /Volumes/prj/Niels_Gehring/newData_March_2018/workflow/mapping/*_STARmapping/ReadsPerGene.out.tab | grep -v 'mate' | cut -f1,4,8,12,16,20,24,28,32,36 | tail -n +5 > newData_March_2018.counts
library('ggplot2')
library('knitr')
library('DESeq2')
library("edgeR")
library('dplyr')
library('RColorBrewer')
library('pheatmap')
library("regionReport")
library('gProfileR')
library("biomaRt")
dir.create('Desktop/Gehring-Casc3/edgeReport',
showWarnings = FALSE, recursive = TRUE)

x <- read.delim("Desktop/Gehring-Casc3/newData_March_2018.counts",
header = F, row.names = 1)

# IH 84667 84669 84671
# WT 84655 84657 84659
# EIF4A3 84661 84663 84665
colnames(x) <- c("IHa", "IHb", "IHc", "WTa", 'WTb', 'WTc', 'EIF4A3a',
'EIF4A3b', 'EIF4A3c')
group <- factor(c(2, 2, 2, 1, 1, 1, 3, 3, 3))
dge <- DGEList(counts = x, group = group)
dge <- calcNormFactors(dge)
plotMDS(dge, col = c(rep("red", 3), rep("blue", 3), rep('green', 3)))

# remove genes with low CPM
# keep <-rowSums(cpm(dge)>=1) >=2
# dge < -dge[keep,]
# d <- cpm(dge)
# # hist(d)
# write.table(dge, "~/Desktop/CPM.txt", sep="\t")

design <- model.matrix(~ 0 + group)
dge <- estimateGLMCommonDisp(dge, design)
# dge <- estimateGLMTrendedDisp(dge, design)
# dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
lrt.IHvsWT <- glmLRT(fit, contrast = c(- 1, 1, 0))
lrt.EIF4A3cvsWT <- glmLRT(fit, contrast = c(0, 1, 1))

sig_IHvsWT <- topTags(lrt.IHvsWT, n = 5000, p.value = 0.05, sort.by = "logFC",
adjust.method = 'fdr')

write.table(sig_IHvsWT,
"Desktop/Gehring-Casc3/top_de_IHvsWT.txt", sep = "\t")

# lrt.IHvsWT$table$FDR <- p.adjust(lrt.IHvsWT$table$PValue,"fdr")

hist(lrt.IHvsWT$table[, 'FDR'], breaks = 20)
hist(lrt.IHvsWT$table[, 'PValue'], breaks = 20)
plotSmear(lrt.IHvsWT,
    de.tags = rownames(
    lrt.IHvsWT$table)[which(lrt.IHvsWT$table$FDR < 0.01)])


edgeReport(
dge, lrt.IHvsWT, project = 'WT vs IH vs EIF4A3', intgroup = 'group',
browse = TRUE, outdir = 'Desktop/Gehring-Casc3/')


# GO term analysis
sel <- row.names(sig_IHvsWT)
universe <- rowSums(cpm(dge) >= 1) >= 2
universe <- as.vector(names(universe[universe]))

goterm <- gprofiler(query = sel, organism = "hsapiens", sort_by_structure = T,
ordered_query = T, significant = T, exclude_iea = T, underrep = F, evcodes = F,
region_query = F, max_p_value = 1, min_set_size = 3, max_set_size = 500, min_isect_size = 0, correction_method = "fdr", hier_filtering = "none", domain_size = "annotated", custom_bg = universe, numeric_ns = "", include_graph = F, src_filter = NULL)
write.table(goterm,
"Desktop/Gehring-Casc3/goterms_lrt.IHvsWT.txt", sep = "\t")

# Create expression file for Enrichment map
normalized_expression_RNAseq <- cpm(dge, normalized.lib.size = TRUE)
keep <- rowSums(cpm(normalized_expression_RNAseq) >= 1) >= 2
normalized_expression_RNAseq <- normalized_expression_RNAseq[keep,]

genenames <- unlist(lapply(rownames(normalized_expression_RNAseq),
function(data) {unlist(strsplit(data, "\\|"))[1]}))

EM_expressionFile_RNAseq <- data.frame(Name = genenames,
normalized_expression_RNAseq)
rownames(EM_expressionFile_RNAseq) <- rownames(normalized_expression_RNAseq)
colnames(EM_expressionFile_RNAseq) <- substring(
colnames(EM_expressionFile_RNAseq), 1, 12)

mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL")
mart = useDataset(mart, dataset = "hsapiens_gene_ensembl")

genes = getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
filters = 'ensembl_gene_id', values = genenames, mart = mart);
genes$description = gsub("\\[Source.*", "", genes$description);


EM_expressionFile_RNAseq <- merge(genes, EM_expressionFile_RNAseq,
all.y = TRUE, by.x = 1, by.y = 1)
colnames(EM_expressionFile_RNAseq)[1] <- "Name"
colnames(EM_expressionFile_RNAseq)[2] <- "Description"
write.table(
    EM_expressionFile_RNAseq,
    'Desktop/Gehring-Casc3/enrichment_map.txt',
    col.name = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)


