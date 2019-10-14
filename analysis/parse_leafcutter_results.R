#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



get.names <- function(x) {
  sec.last.item(str_split(x, '/')[[1]])
}

# load the cluster significance files
cluster_sig_file <-
  Sys.glob('leafcutter/*/*_cluster_significance.txt')
cluster_sig <- lapply(cluster_sig_file, read_tsv)
names(cluster_sig) <- str_split(string = cluster_sig_file,
                                pattern = '/',
                                simplify = TRUE)[, 2]
cluster_sig <-  bind_rows(cluster_sig, .id = 'contrast')

# load the effect sizes files
effec_size_files <- Sys.glob('leafcutter/*/*effect_sizes.txt')
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
write_csv(df, 'leafcutter_junctions.csv')
