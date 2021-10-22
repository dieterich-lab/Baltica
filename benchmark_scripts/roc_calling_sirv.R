library(ROCR)
library(cowplot)

source("load_sirv_data.R")

# df2 <- df
# 
# df2$coordinate <- coordinates 
# df2$comparison <- comparison
# 
# df2 %>% 
#   filter(!is.na(orthogonal)) %>% 
#   distinct(comparison, coordinate) %>% 
#   nrow()
# 
# df2 %>% 
#   filter(!is.na(orthogonal)) %>% 
#   mutate(across(1:4, ~replace_na(.x, 0))) %>% 
#   select(coordinate, comparison, everything()) %>% 
#   write_csv('/prj/Niels_Gehring/baltica_benchmark/1_1_2/input.csv')

length(unique(coordinates[!is.na(df$orthogonal) ]))

compute_prediction <- function(col, ref) {
  col <- df[, col]
  ref <- df[, ref]
  col <- col[!is.na(ref)]
  ref <- ref[!is.na(ref)]
  col <- replace(col, is.na(col), 0)
  prediction(col, ref)
}

pars <- tibble(
  col = c("majiq", "leafcutter", "rmats", "junctionseq"),
  ref = rep("orthogonal", 4) 
) 
prediction <- purrr::pmap(pars, compute_prediction)
names(prediction) <- pars$col
tpr_fpr <- lapply(prediction, performance, measure = "tpr", x.measure = "fpr")
auc <- lapply(prediction, performance, measure = "auc")

method <- unlist(lapply(tpr_fpr, function(x) length(slot(x, "x.values")[[1]])))
auc <- lapply(auc, function(x) round(slot(x, "y.values")[[1]][[1]], 2))

roc_data <- tibble(
  FPR = unlist(lapply(tpr_fpr, function(x) slot(x, "x.values")[[1]])),
  TPR = unlist(lapply(tpr_fpr, function(x) slot(x, "y.values")[[1]])),
  Method = rep(names(method), method),
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Method)) +
  geom_line() +
  geom_point(size=0.5) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme_cowplot(14) +
  theme(legend.position="bottom") +
  scale_color_manual(
    name = NULL,
    labels = setNames(nm=names(auc), paste(names(auc), ' AUCROC=', auc, sep="")),
    values = color_list$method[1:4]) + 
  guides(color=guide_legend(nrow=2,byrow=TRUE))

p

ggsave(
  "../sirv_benchmark/results/roc_calling.pdf", width = 15, height = 15, units = 'cm')
