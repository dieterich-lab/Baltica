source("load_sirv_data.R")

cm <- lapply(
  seq_along(df),
  function(x) {
    
    .x <- df[, c(x, 5)]
    .x <- .x[rowSums(is.na(.x)) !=2, ]
    .x <- replace(.x, is.na(.x), 0)
    
    caret::confusionMatrix(
      factor(.x[[1]] > 0.95),
      factor(.x[[2]] > 0.95),
      positive = "TRUE"
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
