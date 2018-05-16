suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('knitr'))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('RColorBrewer'))
suppressPackageStartupMessages(library('pheatmap'))
suppressPackageStartupMessages(library("regionReport"))
suppressPackageStartupMessages(library('gProfileR'))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("stringr"))

df <- read.csv("DEanalysis/gene_counts.csv",
header = T, row.names = 1)

if (all(startsWith(rownames(df), 'ENSMUS'))) {
    org <- 'mmusculus'
} else if (all(startsWith(rownames(df), 'ENSG'))) {
    org <- 'hsapiens'
}
# else break


group <- factor(str_split(colnames(df), "_", simplify = T)[, 1])
dge <- DGEList(counts = df, group = group)
dge <- calcNormFactors(dge)
plotMDS(dge)

# remove genes with low CPM
keep <- rowSums(cpm(dge) >= 1) >= 2
dge <- dge[keep,]

write.table(dge, "DEanalysis/DGE_kept.txt", sep = "\t")

design <- model.matrix(~ 0 + group)
dge <- estimateGLMCommonDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast = c(- 1, 1))

sig <- topTags(lrt, n = 5000, p.value = 0.05, sort.by = "logFC", adjust.method = 'fdr')
write.table(sig, "DEanalysis/edgeReport/sig_genes.txt", sep = "\t")

plotSmear(lrt, de.tags = rownames(
lrt$table)[which(lrt$table$FDR < 0.01)])

edgeReport(dge, lrt, project = paste(unique(group), sep = " vs "),
intgroup = 'group', browse = TRUE, outdir = 'DEanalysis/edgeReport/')

sel <- row.names(sig)
universe <- rowSums(cpm(dge) >= 1) >= 2
universe <- as.vector(names(universe[universe]))

goterm <- gprofiler(query = sel, organism = org, sort_by_structure = T, ordered_query = T, significant = T, exclude_iea = T, underrep = F, evcodes = F, region_query = F, max_p_value = 1, min_set_size = 3, max_set_size = 500, min_isect_size = 0, correction_method = "fdr", hier_filtering = "none", domain_size = "annotated", custom_bg = universe, numeric_ns = "", include_graph = F, src_filter = NULL)
write.table(goterm, "DEanalysis/goterm.txt", sep = "\t")

# # Create expression file for Enrichment map
# normalized_expression_RNAseq <- cpm(dge, normalized.lib.size = TRUE)
# keep <- rowSums(cpm(normalized_expression_RNAseq) >= 1) >= 2
# normalized_expression_RNAseq <- normalized_expression_RNAseq[keep,]
#
# genenames <- unlist(lapply(rownames(normalized_expression_RNAseq),
# function(data) {unlist(strsplit(data, "\\|"))[1]}))
#
# EM_expressionFile_RNAseq <- data.frame(Name = genenames,
# normalized_expression_RNAseq)
# rownames(EM_expressionFile_RNAseq) <- rownames(normalized_expression_RNAseq)
# colnames(EM_expressionFile_RNAseq) <- substring(
# colnames(EM_expressionFile_RNAseq), 1, 12)
#
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL")
# mart = useDataset(mart, dataset = "hsapiens_gene_ensembl")
#
# genes = getBM(
# attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
# filters = 'ensembl_gene_id', values = genenames, mart = mart);
# genes$description = gsub("\\[Source.*", "", genes$description);
#
#
# EM_expressionFile_RNAseq <- merge(genes, EM_expressionFile_RNAseq,
# all.y = TRUE, by.x = 1, by.y = 1)
# colnames(EM_expressionFile_RNAseq)[1] <- "Name"
# colnames(EM_expressionFile_RNAseq)[2] <- "Description"
# write.table(
# EM_expressionFile_RNAseq,
# 'Desktop/Gehring-Casc3/enrichment_map.txt',
# col.name = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
#
#
