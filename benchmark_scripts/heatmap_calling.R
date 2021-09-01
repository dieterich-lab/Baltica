library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)

df <- readr::read_csv(
  "~/sirv_benchmark/results/SJ_annotated.csv",
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

comparisons <- c("mix3-vs-mix1", "mix2-vs-mix1", "mix3-vs-mix2")
methods <- c("rmats", "junctionseq", "majiq", "leafcutter")

group = str_extract(colnames(mat), str_c(comparisons, collapse = "|"))
method = str_extract(colnames(mat), str_c(methods, collapse = "|"))

mat <- df %>%
  dplyr::select(matches("-vs-")) %>%
  as.matrix(.)

mat[mat > 0.95] <- 1
mat[mat < 0.95] <- 0
mat[is.na(mat)] <- 0

method_cols <- setNames(brewer.pal(n = length(methods), name = "Set1"), methods)
group_cols <- setNames(brewer.pal(n = length(comparisons), name = "Set2"), comparisons)

lengend_method = Legend(
  labels = unique(method), 
  title = "method",
  legend_gp = gpar(fill =  method_cols[unique(method)]))

leg_group = Legend(
  labels = unique(comparisons), 
  title = "group",
  legend_gp = gpar(fill =  group_cols[unique(comparisons)]))


# ra = HeatmapAnnotation(
#   changing = changing,
#   col=list(changing=c('TRUE'='black', 'FALSE'='white')),
#   which = 'row')

col_fun = colorRamp2(c(0, 0.95, 1), c("blue", "white", "red"))

mat <- df %>%
  dplyr::select(matches("-vs-")) %>%
  as.matrix(.) 

mat[is.na(mat)] <- 0

method_cols <- setNames(brewer.pal(n = length(methods), name = "Set1"), methods)
group_cols <- setNames(brewer.pal(n = length(comparisons), name = "Set2"), comparisons)

ca=HeatmapAnnotation(
  method_names=anno_block(
    gp = gpar(fill = method_cols[method]),
    labels = method,
    labels_gp = gpar(col = "white", fontsize = 9)),
  group_names=anno_block(
    gp = gpar(fill = group_cols[group]),
    labels = group,
    labels_gp = gpar(col = "white", fontsize = 6))
)

lengend_method = Legend(
  labels = unique(method), 
  title = "method",
  legend_gp = gpar(fill =  method_cols[unique(method)]))

leg_group = Legend(
  labels = unique(comparisons), 
  title = "group",
  legend_gp = gpar(fill =  group_cols[unique(comparisons)]))

is_SIRV <- startsWith(coordinates, 'SIRV')

ra = HeatmapAnnotation(
  SIRV = is_SIRV,
  col=list(SIRV=c('TRUE'='black', 'FALSE'='white')),
  # changing = changing,
  which = 'row')

col_fun = colorRamp2(c(0, 0.95, 1), c("blue", "white", "red"))

hm=Heatmap(
  mat,
  name='score',
  col = col_fun,
  left_annotation = ra,
  top_annotation = ca,
  cluster_columns = T,
  cluster_rows = F,
  column_split = 1:12,
  # row_split = changing,
  border = TRUE,
  column_title = NULL,
  use_raster=T,
  show_column_names=F)

png("~/sirv_benchmark/heatmap_calling.pdf", width = 8, height = 4, units = "in", res = 200)
draw(hm,
     heatmap_legend_list = list(lengend_method, leg_group))
dev.off()

