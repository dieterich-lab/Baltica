library(tidyverse)

df <- readr::read_csv(
  "results/SJ_annotated.csv",
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
  dplyr::select(matches("-vs-"), "coordinates", "gene_name")

coordinates <- df$coordinates
gene_name <- df$gene_name
df$coordinates <- NULL
df$gene_name <- NULL

## Upset
library(ComplexHeatmap)
library(RColorBrewer)

comparisons <- c("mix3-vs-mix1", "mix2-vs-mix1", "mix3-vs-mix2")
methods <- c("rmats", "junctionseq", "majiq", "leafcutter", "orthogonal")

mat <- df %>%
  dplyr::select(matches("-vs-")) %>%
  as.matrix(.)

mat[mat > 0.95] <- 1
mat[mat < 0.95] <- 0
mat[is.na(mat)] <- 0
mat <- mat[rowSums(mat)>=1,]

comb_mat <- make_comb_mat(mat)

group = str_extract(colnames(mat), str_c(comparisons, collapse = "|"))
method = str_extract(colnames(mat), str_c(methods, collapse = "|"))

comb_mat <- comb_mat[comb_degree(comb_mat) > 1]

p2 <- UpSet(
  comb_mat,
  comb_order = order(rev(comb_degree(comb_mat))),
  pt_size = unit(2, "mm"), lwd = 1.5,
  right_annotation = rowAnnotation(
    group = group,
    method = method,
    col = list(
      group = setNames(brewer.pal(n = length(comparisons), name = "Set2"), comparisons) ,
      method = setNames(brewer.pal(n = length(methods), name = "Set1"), methods)
    )
  )
)


p2 <- UpSet(
  comb_mat,
  comb_order = order(rev(comb_degree(comb_mat))),
  pt_size = unit(2, "mm"), lwd = 1.5,
  right_annotation = rowAnnotation(
    group = group,
    method = method,
    col = list(
      group = setNames(brewer.pal(n = length(comparisons), name = "Set2"), comparisons) ,
      method = setNames(brewer.pal(n = length(methods), name = "Set1"), methods)
    )
  )
)
pdf(
   "upset_calling.pdf",
   width = 10,
   height = 5,
   useDingbats = FALSE
 )
p2
dev.off()

