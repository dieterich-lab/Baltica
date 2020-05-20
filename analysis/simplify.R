suppressPackageStartupMessages({
  library(GenomicRanges)
  library(stringr)
  library(igraph)
})
MISSING file_names
files <- stringr::str_split(file_names, pattern = ',', simplify = TRUE)

df <- apply(files, 2, read.csv, stringsAsFactors = FALSE)

names(df) <- stringr::str_split(files, pattern='/', simplify = TRUE)[,1]

gr <- lapply(df, function(x) {
    .tmp <- GRanges(x)
    mcols(.tmp)  <- NULL
    .tmp
})

gr <- stack(as(gr, 'GRangesList'))

hits <- findOverlaps(gr, drop.self = TRUE)
query <- gr[queryHits(hits)]
subject <- gr[subjectHits(hits)]
# filter the overlaps by max 2 nt difference in start and stop
start_dif <- abs(start(query) - start(subject))
end_dif <- abs(end(query) - end(subject))
tolerance <- 2
hits <- hits[start_dif <= tolerance & end_dif <= tolerance]
# uses the quick union algorithm to find the isolates
hits <- as.data.frame(hits)
hits_graph <- igraph::graph.data.frame(hits)
membership <- igraph::clusters(hits_graph)$membership
hits <-
  merge(
    hits,
    stack(membership),
    by.x = "queryHits",
    by.y = "ind",
    all.x = TRUE
  )
# use the membership id to group isolates
gr$SJ_membership <- NA
gr[hits$queryHits]$SJ_membership <- hits$values
# uses numbers higher the number of groups to name the junctions groups
singletons <- cumsum(is.na(gr$SJ_membership))
singletons <- singletons + max(gr$SJ_membership, na.rm = T)
gr$SJ_membership[is.na(gr$SJ_membership)] <-
  singletons[is.na(gr$SJ_membership)]

.cols  <- intersect(colnames(df$majiq), colnames(df$leafcutter))

.cols <- .cols[3: 16]

df <- lapply(df, `[`, .cols)

df <- do.call("rbind", df)

as_type <- apply(
    stringr::str_split(df$as_type, ';', simplify = T),
    1,
    function(x){ paste0(setdiff(unique(x), c('JS', 'JE', 'NA', '')), collapse = ';') }
)

df$as_type <- as_type

df$SJ_membership <- gr$SJ_membership

.cols2  <- .cols[!.cols %in% c('method', 'comparison') ]

df_method  <- aggregate(method ~ comparison + SJ_membership,
                        data = df,
                        FUN = function(x) paste0(unique(x), collapse = ';'))

df_cols <- aggregate(as.matrix(df[,.cols2]) ~ comparison + SJ_membership,
                        data = df,
                        FUN = function(x) x[[1]])

df2  <-  merge(df_method, df_cols, by = c('comparison', 'SJ_membership') )

df2 <- df2[order(df2$chr, df2$start), ]

df2 <- subset(df2, select=-c(SJ_membership, is_canonical))