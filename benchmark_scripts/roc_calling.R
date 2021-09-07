library(tidyverse)

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
) %>%
  dplyr::select(matches("-vs-")) %>% 
  # NA are not called
  mutate( 
    across(everything(), ~ replace_na(.x, 0))
  ) %>% # calling task
  mutate_all( 
    list(~ if_else(. > 0.05, T, F))
  ) %>%
  pivot_longer(everything(),
               names_to = c(".value", 'comparison'),
               names_pattern = "(.+)_(.+)"
  ) 

df$comparison <- NULL

compute_cm <- function(col, ref) {
  case_when(
    df[[col]] == 0 & df[[ref]] == 0 ~ "NA",
    df[[col]] > 0.95 & df[[ref]] == 1 ~ "TP",
    df[[col]] > 0.95 & df[[ref]] != 1 ~ "FP",
    df[[col]] <= 0.95 & df[[ref]] != 1 ~ "TN",
    df[[col]] <= 0.95 & df[[ref]] == 1 ~ "FN",
  )
}

pars <- tibble(
  col = c("majiq", "leafcutter", "rmats", "junctionseq"),
  ref = rep("orthogonal", 4) 
) 

cm_result <- purrr::pmap(pars, compute_cm)
names(cm_result) <- pars$col

cm_result %>%
  purrr::map(~ table(.)) %>%
  bind_rows(., .id = "sample")

# status_leafcutter_mix2vsmix1 = compute_cm("leafcutter_mix2-vs-mix1", "orthogonal_result_mix2-vs-mix1")

library(ROCR)
library(ggplot2)
library(cowplot)

compute_prediction <- function(col, ref) {
  .x <- df[, c(col, ref)]
  .x <- filter_all(.x, any_vars(. != 0))
  prediction(.x[[col]], 1 - .x[[ref]])
}
prediction <- purrr::pmap(pars, compute_prediction)
names(prediction) <- pars$col
tpr_fpr <- lapply(prediction, performance, measure = "tpr", x.measure = "fpr")
auc <- lapply(prediction, performance, measure = "auc")

method <- unlist(lapply(tpr_fpr, function(x) length(slot(x, "x.values")[[1]])))
round(unlist(lapply(auc, function(x) slot(x, "y.values")[[1]])), 2)
auc <- lapply(auc, function(x) round(slot(x, "y.values")[[1]][[1]], 2))


roc_data <- tibble(
  FPR = unlist(lapply(tpr_fpr, function(x) slot(x, "x.values")[[1]])),
  TPR = unlist(lapply(tpr_fpr, function(x) slot(x, "y.values")[[1]])),
  Method = rep(paste0(names(method), " (AUC=", auc[names(method)], ")"), method)
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Method)) +
  geom_line() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme_cowplot(16) +
  labs(title = "SIRV Benchmark - Calling") +
  theme(legend.position = c(0.60, 0.15))

ggsave("sirv_benchmark_roc_calling.pdf")
