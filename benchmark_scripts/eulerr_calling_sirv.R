library(eulerr)

source("~/Baltica/benchmark_scripts/load_sirv_data.R")

# keep <- !is.na(df$orthogonal)

mat <- df %>%
  mutate(
    across(everything(), ~ replace_na(.x, 0)),
    comparison = NULL
  ) %>%
  as.matrix(df)

# mat <- mat[keep, ]  
x[x > 0.95] <- 1
x[x < 0.95] <- 0

# gene_name <- gene_name[keep]
# coordinates <- coordinates[keep]

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


input_sj <- lapply(
  as.data.frame(mat),
  function(i) unique(coordinates[as.logical(i)])
)

input_genes <- lapply(
  as.data.frame(mat),
  function(i) unique(na.omit(gene_name[as.logical(i)]))
)

pdf('../sirv_benchmark/results/eulerr_gene_sirv_calling.pdf')
plot(
  euler(input_genes[1:5], shape = "ellipse"),
  fills = color_list$method[names(input_genes[1:5])],
  edges = F,
  fontsize = 14,
  quantities = list(fontsize = 14),
  labels = FALSE,
  legend = list(nrow = 1, ncol = 5, side = "bottom", fontsize = 16)
)
dev.off()

pdf('../sirv_benchmark/results/eulerr_sj_sirv_calling.pdf')
plot(
  euler(input_sj[1:5], shape = "ellipse"),
  fills = color_list$method[names(input_sj[1:5])],
  edges = F,
  # adjust_labels = TRUE,
  quantities = list(fontsize = 10),
  labels = FALSE,
  legend = list(nrow = 1, ncol = 5, side = "bottom", fontsize = 16))

dev.off()
