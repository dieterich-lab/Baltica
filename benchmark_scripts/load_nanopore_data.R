library(tidyverse)

df <- readr::read_csv(
  "../nanopore_benchmark/results/SJ_annotated.csv",
  col_types = cols(
    .default = col_double(),
    is_novel = col_logical(),
    gene_name = col_character(),
    transcript_name = col_character(),
    class_code = col_character(),
    exon_number = col_character(),
    coordinates = col_character()
  )
) 

gene_name <- df$gene_name
coordinates <- df$coordinates
is_novel <- df$is_novel

colnames(df) <- gsub(x=colnames(df), pattern = "_SMG7KO2SMG6KD-vs-control", '')
colnames(df) <- gsub(x=colnames(df), pattern = "_NA", '')

methods <- c("rmats", "junctionseq", "majiq", "leafcutter", "orthogonal")
df <- df %>% dplyr::select(matches(methods)) 

color_list <- list(
  method = setNames(
    scales::brewer_pal(palette = "Set2")(length(methods)),
    methods
  )
)
