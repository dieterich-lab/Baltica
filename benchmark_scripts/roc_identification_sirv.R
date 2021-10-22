source("load_sirv_data.R")

library(ROCR)
library(ggplot2)
library(cowplot)

df <- df %>% 
  filter(!is.na(orthogonal)) %>% 
  mutate(across(1:4, ~ifelse(is.na(.), yes=0, 1)))

compute_prediction <- function(col, ref) {
  .x <- df[, c(col, ref)]
  prediction(.x[[col]], .x[[ref]])
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
  Method = rep(names(method), method)
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

ggsave("../sirv_benchmark/results/roc_identification.pdf", width = 15, height = 15, units = 'cm')
register(MulticoreParam(workers=3))

