source("load_sirv_data.R")

df <- df %>% 
  replace(is.na(.), 0) %>% 
  mutate_all(~ifelse(. < 0.95, 0, 1))

cm <- lapply(
  seq_along(df),
  function(x) {
    caret::confusionMatrix(
      as.factor(df[[x]]),
      as.factor(df$orthogonal),
      positive = "1"
    )
  }
)

names(cm) <- colnames(df)

lapply(
  names(cm),
  function(x) {
    writeLines(
      capture.output(
        cm[[x]]
      ),
      paste0("../sirv_benchmark/results/confmat_", x, ".tab")
    )
  }
)
