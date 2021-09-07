library(tidyverse)
library(reshape2)

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
  dplyr::select(matches("-vs-"), "coordinates", "gene_name")

coordinates <- df$coordinates
gene_name <- df$gene_name
df$coordinates <- NULL
df$gene_name <- NULL
df <- df %>% mutate(
  across(everything(), ~ replace_na(.x, 0))
)

cormat <- apply(df, 2, function(x){ -log10(x+1e-10) })
cormat <- round(
  cor(cormat, method='pearson',  use="pairwise.complete.obs"),
  2)
cormat[lower.tri(cormat)]=NA
cormat=melt(cormat)
cormat=cormat[!is.na(cormat$value),]
ggplot(cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.7, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
ggsave('../sirv_benchmark/results/heatmap_pearson.png')
