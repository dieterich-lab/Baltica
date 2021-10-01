source("load_sirv_data.R")

library(ComplexHeatmap)
library(RColorBrewer)

set.seed(123)

colnames(df) <- gsub(x=colnames(df), "orthogonal", 'SIRV') 
methods <- c("rmats", "junctionseq", "majiq", "leafcutter", "SIRV")
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
new_idx <- order(mat[, 5])
mat <- mat[new_idx, ]

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

png("~/Baltica/sirv_benchmark/results/heatmap_calling.png", width = 8, height = 4, res=300, units = 'in')
draw(hm,
     heatmap_legend_list = list(lengend_method))
dev.off()


