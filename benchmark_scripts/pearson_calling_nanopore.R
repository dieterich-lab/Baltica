library(reshape2)

source("load_nanopore_data.R")

df <- df %>% 
  filter(!is.na(orthogonal)) %>% 
  replace(is.na(.), 0) 

colnames(df) <- gsub(x=colnames(df), "orthogonal", 'nanopore') 
cormat <- apply(df, 2, function(x){ -log10(x+1e-10) })
cormat <- as.matrix(df)
cormat <- round(
  cor(cormat, method='pearson',  use="pair†wise.complete.obs"),
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
    legend.position = c(0.7, 0.78),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
ggsave('../nanopore_benchmark/results/heatmap_pearson.pdf')
