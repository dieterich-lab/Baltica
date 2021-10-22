source("load_nanopore_data.R")

library(ComplexHeatmap)
library(RColorBrewer)

set.seed(123)

colnames(df) <- gsub(x=colnames(df), "orthogonal", 'nanopore') 
methods <- c("rmats", "junctionseq", "majiq", "leafcutter", "nanopore")
methods <- str_extract(colnames(df), str_c(methods, collapse = "|"))
methods <- methods[!is.na(methods)]

df <- df %>% 
  filter(rowSums(., na.rm=TRUE) > 0, ) %>% 
  filter(!is.na(nanopore)) %>% 
  replace(is.na(.), 0) 

mat <- as.matrix(df)

idx <- data.frame(
  is=ifelse(df$nanopore>0.95, 1, 0)
)

idx <- idx %>%
  mutate(n=row_number()) %>% 
  group_by(is) %>% 
  sample_n(500)

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

col_fun = circlize::colorRamp2(c(0, 0.95, 1), c("blue", "white", "red"))

hm=Heatmap(
  mat,
  name='score',
  col = col_fun,
  top_annotation = ca,
  cluster_columns = F,
  cluster_rows = F,
  column_split = 1:5,
  border = TRUE,
  column_title = NULL,
  show_column_names=F)


png(
  "~/Baltica/nanopore_benchmark/results/heatmap_calling.png", 
  width = 8, height = 4, res=300, units = 'in')
draw(hm,
     heatmap_legend_list = list(lengend_method))
dev.off()


