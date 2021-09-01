library(tidyverse)

# input
# output
# method = c("majiq", "leafcutter", "rmats", "junctionseq"),
# group = c("mix2-vs-mix1", "mix3-vs-mix1", "mix3-vs-mix2")

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

compute_cm <- function(col, ref) {
  case_when(
    df[[col]] == 0 & df[[ref]] == 0 ~ "NA",
    df[[col]] > 0.95 & df[[ref]] == 1 ~ "TP",
    df[[col]] > 0.95 & df[[ref]] != 1 ~ "FP",
    df[[col]] <= 0.95 & df[[ref]] != 1 ~ "TN",
    df[[col]] <= 0.95 & df[[ref]] == 1 ~ "FN",
  )
}

pars <- expand_grid(
  method = c("majiq", "leafcutter", "rmats", "junctionseq"),
  group = c("mix2-vs-mix1", "mix3-vs-mix1", "mix3-vs-mix2")
)
pars <- pars %>%
  mutate(col = str_glue("{method}_{group}"))
pars <- pars %>%
  mutate(ref = str_glue("orthogonal_{group}"))
pars$method <- NULL
pars$group <- NULL

cm_result <- purrr::pmap(pars, compute_cm)
names(cm_result) <- pars$col
table(cm_result$`majiq_mix2-vs-mix1`)

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
