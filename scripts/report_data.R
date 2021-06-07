#! /usr/bin/env Rscript
## ---------------------------
## Script name: report_data
## Purpose of script: Provide data to baltica report
## Author: Thiago Britto-Borges
## Date Created: 2021-05-20
## Copyright (c) Thiago Britto-Borges, DieterichLab 2021
## Email: thiago.brittoborges@uni-heidelberg.de
suppressPackageStartupMessages({
    library(tidyverse)
    library(GenomicRanges)
    library(yaml)
    library(rtracklayer)
    library(ggplot2)
})

config <- yaml::read_yaml("config.yaml")
gtf <- rtracklayer::import(config$ref)
genes <- subset(gtf, type == "gene" & gene_biotype == "protein_coding")

sj_out_files <- lapply(
    config$samples,
    function(x) {
        x <- file.path(config$sample_path, x)
        x <- sub(x = x, "Aligned.noS.bam", "SJ.out.tab")
    }
)

sj_out <- lapply(
    sj_out_files,
    read_delim,
    delim = "\t",
    col_names = c(
        "chr", "start", "end", "strand", "intron_motif", "annotation",
        "uniq_map_read_count", "multi_map_read_count", "max_overhang"
    ),
    col_types = cols(chr = col_character(), .default = col_integer())
)

sj_out <- bind_rows(sj_out, .id = "sample")
sj_out$strand <- case_when(
    sj_out$strand == 0 ~ "*",
    sj_out$strand == 1 ~ "+",
    sj_out$strand == 2 ~ "-"
)

# : 0: noncanonical, 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
sj_out$intron_motif <- case_when(
    sj_out$intron_motif == 0 ~ "noncanonical",
    sj_out$intron_motif == 1 ~ "GT/AG",
    sj_out$intron_motif == 2 ~ "CT/AC",
    sj_out$intron_motif == 3 ~ "GC/AG",
    sj_out$intron_motif == 4 ~ "CT/GC",
    sj_out$intron_motif == 5 ~ "AT/AC",
    sj_out$intron_motif == 6 ~ "GT/AT"
)

# 0: unannotated, 1: annotated, only if an input gene annotations file was used
sj_out$annotation <- case_when(
    sj_out$annotation == 0 ~ "unannotated",
    sj_out$annotation == 1 ~ "annotated",
)


sj_out <- sj_out %>%
    mutate(range = str_glue_data(., "{chr}:{start}-{end}:{strand}"))

plots_hist_total <- sj_out %>%
    ggplot(aes(x = uniq_map_read_count + multi_map_read_count)) +
    facet_wrap(~sample, ncol = 4) +
    geom_histogram(alpha = .3) +
    scale_x_log10() +
    xlab("log10(total junction reads)") +
    ylab("freq") +
    theme_bw(16)

## PCA ####
pca <- sj_out %>%
    pivot_wider(
        id_cols = c(chr, start, end, strand),
        values_from = uniq_map_read_count,
        names_from = sample,
        values_fill = 0
    ) %>%
    select(-c(chr, start, end, strand))

groups <- colnames(pca)
pca <- as.matrix(pca)
rv <- matrixStats::rowVars(pca)
top <- order(rv, decreasing = TRUE)[seq_len(1000)]
pca <- pca[top, ]

# mode(pca) = "numeric"
pca <- t(pca)
pca <- prcomp(pca, scale. = TRUE, center = TRUE)

percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- paste0(
    colnames(pca$x),
    "(", paste(as.character(percentage), "%", ")",
        sep =
            ""
    )
)
pca <- as.data.frame(pca$x)
pca$group <- str_split(groups, "_", simplify = T)[, 1]

plot_pca <- ggplot(pca, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 10) +
    theme_minimal() +
    xlab(percentage[1]) +
    ylab(percentage[2]) +
    theme_bw(16) +
    labs(title = "Junction reads", subtitle = "Top 1e4 variance junctions")


res <- read_csv(
    "results/SJ_annotated_assigned.csv",
    col_types = cols(
        seqnames = col_character(),
        start = col_double(),
        end = col_double(),
        width = col_double(),
        strand = col_character(),
        is_novel = col_character(),
        gene_name = col_character(),
        transcript_name = col_character(),
        class_code = col_character(),
        exon_number = col_character(),
        comparison = col_character(),
        method = col_character(),
        score = col_character(),
        as_type = col_character()
    )
)


url <- str_glue(
    "http://genome.ucsc.edu/cgi-bin/hgTracks?db={config$assembly}hg38&position="
)
res <- res %>%
    mutate(coordinate = str_glue_data(., "{seqnames}:{start}-{end}:{strand}")) %>%
    mutate(coordinate = paste0("<a href='", url, str_sub(coordinate, 1, -3), "' target='_blank'>", coordinate, "</a>")) %>%
    select(-c(seqnames, start, end, width, strand, exon_number)) %>%
    select(c(coordinate, comparison, method, score, gene_name, transcript_name, class_code, as_type, is_novel))

simplify <- function(x, remove = c()) {
    x %>%
        str_split(";") %>%
        lapply(., function(x) {
            setdiff(x, remove)
        }) %>%
        lapply(., paste0, collapse = ";") %>%
        unlist()
}

for (col in c("is_novel", "gene_name", "transcript_name", "class_code", "comparison", "method", "as_type")) {
    if (col == "as_type") {
        res[[col]] <- simplify(res[[col]], remove = c("JS", "JE"))
        next
    }

    res[[col]] <- simplify(res[[col]])
}

if (file.exists("leafcutter_junctions_complete.csv")) {
    leafcutter_complete <- read_csv("leafcutter_junctions_complete.csv")
} else {
    system("parse_leafcutter_output.R -c 1.1 -o leafcutter_junctions_complete.csv",
        wait = T
    )
    leafcutter_complete <- read_csv("leafcutter_junctions_complete.csv")
}

leafcutter_complete <- leafcutter_complete %>%
    mutate(is_sig = abs(deltapsi) > 0.1 & p.adjust < 0.05)

gr_leafcutter_complete <- leafcutter_complete %>%
    mutate_at(vars(chr), replace_na, "0") %>%
    makeGRangesFromDataFrame(.)

hits <- findOverlaps(genes, gr_leafcutter_complete)
hits <- as.data.frame(hits)
hits$gene_name <- genes[hits$queryHits, ]$gene_name
hits <- aggregate(gene_name ~ subjectHits, data = hits, c)
hits$gene_name <- lapply(hits$gene_name, unique)
hits$gene_name <- lapply(hits$gene_name, paste0, collapse = ";")
hits$gene_name <- unlist(hits$gene_name)
leafcutter_complete[hits$subjectHits, ]$genes <- unlist(hits$gene_name)

top_rows <- leafcutter_complete %>%
    rownames_to_column() %>%
    group_by(comparison) %>%
    filter(is_sig) %>%
    top_n(10, abs(deltapsi)) %>%
    pull(rowname) %>%
    as.numeric(.)

save.image("baltica_report.RData")