source("load_nanopore_data.R")

cm <- lapply(
  seq_along(df),
  function(x) {
    caret::confusionMatrix(
      as.factor(ifelse(df[[x]] > 0.95, 1, 0)),
      as.factor(ifelse(df$orthogonal > 0.95, 1, 0)),
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
      paste0("../nanopore_benchmark/results/confmat_", x, ".tab")
    )
  }
)