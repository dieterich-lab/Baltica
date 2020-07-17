#!/usr/bin/env Rscript
# Title     : annotate_SJ from many DJU methods
# Objective : To annotate SJ called DS by many DJU methods with gene and transcript level information
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
  library(openxlsx)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result",,
    metavar = "character",
    default = "results/SJ_annotated_assigned.xlsx"
  ),
    make_option(
    c("-o", "--output"),
    type = "character",
    help = "Path to the output file",
    metavar = "character",
    default = "results/SJ_annotated_assigned_simple.xlsx"

  )
)

opt <- parse_args(OptionParser(option_list = option_list))
x <- read_csv(opt$input)

simplify <- function(x, remove=c()){   
    x  %>% str_split(';') %>%
    lapply(., function(x) { setdiff(x, remove) }) %>% 
    lapply(., paste0, collapse = ';')   %>% 
    unlist()
}

x <- x %>% 
    select(-c(X1, Row.names, X, width)) %>% 
    mutate(range = as.character(str_glue('{seqnames}:{start}-{end}:{strand}'))) %>% 
    select(-c(seqnames, start, end, strand))

for (col in c('gene_name', 'transcript_name', 'class_code', 'exon_number', 'comparison', 'method', 'as_type' )){
    if (col == 'as_type'){
        x[[col]]  <- simplify(x[[col]], remove = c('JS', 'JE'))    
        next
    }
    
    x[[col]]  <- simplify(x[[col]])
    
}


openxlsx::write.xlsx(
  x,
  file = opt$output,
  sheetName = "Sheet1",
)
