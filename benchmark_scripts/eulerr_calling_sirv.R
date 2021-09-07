library(eulerr)

source("~/Baltica/benchmark_scripts/load_sirv_data.R")

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

fit <- euler(mat, shape = "ellipse")

plot(
  fit,
  fills = color_list$method,
  edges = F,
  fontsize = 14,
  quantities = list(fontsize = 16),
  labels = FALSE,
  legend = list(nrow = 1, ncol = 5, side = "bottom", fontsize = 16),
  main = "Intesection of SJ per method"
)

library(VennDiagram)

input_sj <- lapply(
  as.data.frame(mat),
  function(i) unique(coordinates[i])
)

venn.diagram(
  x = as.list(mat),
  category.names = colnames(mat),
  filename = "venn_sirv_calling.png"
  #   output=TRUE
)

input <- lapply(as.data.frame(mat), function(i) unique(gene_name[i]))

input$junctionseq

venn.diagram(
  x = input,
  #   category.names = colnames(mat),
  filename = "venn_sirv_calling.png",
  output = TRUE
)
