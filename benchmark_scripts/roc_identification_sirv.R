source("load_sirv_data.R")

library(ROCR)
library(ggplot2)
library(cowplot)

# NA are not called
df <- df %>% 
  mutate( 
    across(everything(), ~ifelse(is.na(.), 0, 1)),
    comparison = NULL
  ) 

compute_prediction <- function(col, ref) {
  .x <- df[, c(col, ref)]
  .x <- filter_all(.x, any_vars(. != 0))
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
  Method = rep(paste0(names(method), " (AUC=", auc[names(method)], ")"), method)
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Method)) +
  geom_line() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme_cowplot(14) +
  labs(title = "SIRV benchmark - identification task") +
  theme(legend.position = c(0.60, 0.3))
p

ggsave("../sirv_benchmark/results/roc_identification.pdf")

