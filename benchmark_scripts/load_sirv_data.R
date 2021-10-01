library(tidyverse)

df <- readr::read_csv(
  "../sirv_benchmark/results/SJ_annotated.csv",
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

df <- df %>% 
  # pivot to remove the mixes
  pivot_longer(matches("-vs-"),
               names_to = c(".value", 'comparison'),
               names_pattern = "(.+)_(.+)"
  ) 

gene_name <- df$gene_name
coordinates <- df$coordinates
# is_novel <- df$is_novel

methods <- c("rmats", "junctionseq", "majiq", "leafcutter", "orthogonal")

df <- df %>% dplyr::select(matches(methods)) 

color_list <- list(
  method = setNames(
    scales::brewer_pal(palette = "Set2")(length(methods)),
    methods
  )
)
