source("load_sirv_data.R")

library(ComplexHeatmap)
library(RColorBrewer)

set.seed(123)

methods <- c("rmats", "junctionseq", "majiq", "leafcutter", "orthogonal")
methods <- str_extract(colnames(df), str_c(methods, collapse = "|"))
methods <- methods[!is.na(methods)]

mat <- df %>% 
  mutate( 
    across(everything(), ~ replace_na(.x, 0)),
    comparison = NULL
  ) %>% 
  as.matrix(df)

is_SIRV <- startsWith(coordinates, 'SIRV')

idx <- data_frame(
  is_SIRV=is_SIRV)

idx <- idx %>%
  mutate(n=row_number()) %>% 
  group_by(is_SIRV) %>% 
  sample_n(504)

method <- str_extract(colnames(mat), str_c(methods, collapse = "|"))
method_cols <- setNames(brewer.pal(n = length(methods), name = "Set2"), methods)

mat <- mat[idx$n, ]
ca=HeatmapAnnotation(
  method_names=anno_block(
    gp = gpar(fill = method_cols[method]),
    labels = method,
    labels_gp = gpar(col = "white", fontsize = 9))
)

lengend_method = Legend(
  labels = unique(method),
  title = "method",
  legend_gp = gpar(
    fontsize = 14,
    fontface = "bold",
    fill =  method_cols[unique(method)]))

ra = HeatmapAnnotation(
  SIRV = is_SIRV[idx$n],
  col=list(SIRV=c('TRUE'='black', 'FALSE'='white')),
  which = 'row')

col_fun = circlize::colorRamp2(c(0, 0.95, 1), c("blue", "white", "red"))

hm=Heatmap(
  mat,
  name='score',
  col = col_fun,
  left_annotation = ra,
  top_annotation = ca,
  cluster_columns = F,
  cluster_rows = F,
  column_split = 1:5,
  row_split = is_SIRV[idx$n],
  border = TRUE,
  column_title = NULL,
  show_column_names=F)

png("~/Baltica/sirv_benchmark/results/heatmap_calling.pdf", width = 8, height = 4, units = "in", res = 200)
draw(hm,
     heatmap_legend_list = list(lengend_method))
dev.off()


