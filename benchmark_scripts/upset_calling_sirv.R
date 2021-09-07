source("load_sirv_data.R")

library(ComplexHeatmap)
library(RColorBrewer)

methods <- c("rmats", "junctionseq", "majiq", "leafcutter", "orthogonal")
methods <- str_extract(colnames(df), str_c(methods, collapse = "|"))
methods <- methods[!is.na(methods)]

mat <- df %>% 
  mutate( 
    across(everything(), ~ replace_na(.x, 0)),
    comparison = NULL
  ) %>% 
  as.matrix(df)

mat[mat > 0.95] <- 1
mat[mat < 0.95] <- 0
mat[is.na(mat)] <- 0
mat <- mat[rowSums(mat) >= 1, ]
comb_mat <- make_comb_mat(mat)
method <- str_extract(colnames(mat), str_c(methods, collapse = "|"))
color_list <- list(
  method = setNames(
    scales::brewer_pal(palette = "Set2")(length(methods)),
    methods
  )
)
comb_mat <- comb_mat[comb_degree(comb_mat) > 1]

UpSet(
  comb_mat,
  comb_order = order(rev(comb_degree(comb_mat))),
  top_annotation=upset_top_annotation(comb_mat, add_numbers = T),
  right_annotation = upset_right_annotation(
    comb_mat, add_numbers = T),
  lwd = 0.8,
  pt_size = unit(2, "mm"))

pdf(
  "../sirv_benchmark/results/upset_calling_sirv.pdf",
  width = 6,
  height = 4,
  useDingbats = FALSE
)
draw(us)
dev.off()
