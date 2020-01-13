#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(tidylog)
  library(optparse)
})


option_list = list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "majiq/voila/*.tsv",
    help = "Path to Majiq result files",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "majiq/majiq_junctions.csv",
    help = "Output file",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

files <- Sys.glob(opt$input)

# rename column names from majiq result due to the presence of spaces
col_names <- "Gene_Name
Gene_ID
LSV_ID
E_dPSI_per_LSV_junction
P_dPSI_beq_per_LSV_junction
P_dPSI_leq_per_LSV_junction
ref_E_PSI
alt_E_PSI
LSV_Type
A5SS
A3SS
ES
Num_Junctions
Num_Exons
chr
strand
Junctions_coords
Exons_coords
IR_coords
UCSC_LSV_Link
De_Novo_Junctions"

read_majiq_out <- function(x) {
  read_tsv(
    x,
    col_names = strsplit(col_names, '\n')[[1]],
    comment = '#',
    cols(.default = col_character(),
         A5SS = col_logical(),
         A3SS = col_logical(),
         ES = col_logical(),
         Num_Junctions = col_double(),
         Num_Exons = col_double()
         )
    )
}

res <- lapply(files, read_majiq_out)
names(res) <- lapply(files, function(x){ tools::file_path_sans_ext(basename(x)) })

res <- bind_rows(res, .id='comparison') %>%
  separate_rows(
    "Junctions_coords",
    "E_dPSI_per_LSV_junction",
    "P_dPSI_beq_per_LSV_junction",
    "P_dPSI_leq_per_LSV_junction",
    "ref_E_PSI",
    "alt_E_PSI",
    sep=';',
    convert = T)

write_csv(res, opt$output)
