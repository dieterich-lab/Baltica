#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
    })

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) > 1) {
  stop("No more than 1 argument should be supplied.n", call.=FALSE)
} else if (!dir.exists(args[1])){
  stop("Input directory not found.n", call.=FALSE)
}

out.path <- args[1]

# load the cluster significance files
cluster_sig_file <- Sys.glob(
  file.path(
    out.path,
    '/*/*_cluster_significance.txt')
  )
cluster_sig <- lapply(cluster_sig_file, read_tsv)
names(cluster_sig) <- str_split(string = cluster_sig_file,
                                pattern = '/',
                                simplify = TRUE)[, 2]
cluster_sig <-  bind_rows(cluster_sig, .id = 'contrast')

# load the effect sizes files
effec_size_files <- Sys.glob(
  file.path(
    out.path,
    '/*/*effect_sizes.txt')
  )
es <- lapply(effec_size_files, read_tsv)
names(es) <- str_split(effec_size_files, '/', simplify = TRUE)[,  2]
es <- bind_rows(es, .id = 'contrast')

# parse the intron column for merging
es$intron <- str_replace_all(es$intron, '_', ':')

intron <- read_delim(
  es$intron,
  delim = ':',
  col_names = c('chr', 'start', 'end', 'clu', 'clu_number', 'strand')
)
es <- bind_cols(es, intron)
# cluster will be the pivot for merging
es$cluster <- as.character(str_glue_data(es, "{chr}:{clu}_{clu_number}_{strand}"))
# merge the files
df <- inner_join(es, cluster_sig, by = c('contrast', 'cluster'))
df$chr <- gsub('chr', '', df$chr)
df <- select(df, -c('clu', 'clu_number'))
# create a unique junction column for each row
message('Writting the parsed output to ', file.path(out.path, 'leafcutter_junctions.csv'))
write_csv(df, file.path(out.path, 'leafcutter_junctions.csv'))
