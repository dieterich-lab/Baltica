#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(optparse)
  library(GenomicRanges)
})


args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (Leafcutter output directory).\n",
       call. = FALSE)
} else if (length(args) > 1) {
  stop("No more than 1 argument should be supplied.\n", call. = FALSE)
} else if (!dir.exists(args[1])) {
  stop("Leafcutter output directory not found.\n", call. = FALSE)
}

out.path <- args[1]

message("Loading the cluster significance files")
cluster_sig_file <- Sys.glob(file.path(out.path, '/*/*_cluster_significance.txt'))

cluster_sig <- lapply(
  cluster_sig_file,
  read_tsv,
  col_types = c(
    cluster = col_character(),
    status = col_character(),
    loglr = col_double(),
    df = col_double(),
    p = col_double(),
    p.adjust = col_double(),
    genes = col_character()
  )
)

file_names <- str_split(string = cluster_sig_file,
                        pattern = '/',
                        simplify = TRUE)

names(cluster_sig) <- file_names[, ncol(file_names) - 1]
cluster_sig <- bind_rows(cluster_sig, .id = 'comparison')

message("Loading the effect sizes files")
effec_size_files <- Sys.glob(file.path(out.path,
                                       '/*/*effect_sizes.txt'))
es <- lapply(
  effec_size_files,
  read_tsv,
  col_types = c(
    .default = col_double(),
    intron = col_character(),
    logef = col_double(),
    deltapsi = col_double()
  )
)

es_file_names <- str_split(string = effec_size_files, pattern = '/', simplify = TRUE)
names(es) <- es_file_names[, ncol(es_file_names) - 1]
es <- bind_rows(es, .id = 'comparison')

# parse the intron column for merging
es$intron <- str_replace_all(es$intron, '_', ':')

intron <- read_delim(es$intron, delim = ':', col_names = c('chr', 'start', 'end', 'clu', 'clu_number', 'strand')
)
es <- bind_cols(es, intron)
# cluster will be the pivot for merging
es$cluster <- as.character(str_glue_data(es, "{chr}:{clu}_{clu_number}_{strand}"))
message("Merging tables")

res <- inner_join(es, cluster_sig, by = c('comparison', 'cluster'))
res$chr <- gsub('chr', '', res$chr)
# add ranks for psi per cluster
res <- res %>%
  arrange(contrast, cluster, ref_psi) %>%
  group_by(contrast, cluster) %>%
  mutate(ref_rank = rank(ref_psi, ties.method = "average")) %>%
  ungroup() %>%
  arrange(contrast, cluster, alt_psi) %>%
  group_by(contrast, cluster) %>%
  mutate(alt_rank = rank(alt_psi, ties.method = "average")) %>%
  ungroup()

res <- select(res, -c('clu', 'clu_number'))
# create a unique junction column for each row
message('Writting the parsed output to ', file.path(out.path, 'leafcutter_junctions.csv'))
write_csv(res, file.path(out.path, 'leafcutter_junctions.csv'))
