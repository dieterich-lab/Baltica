#!/usr/bin/env Rscript
suppressPackageStartupMessages({
                                 library(tidyr)
                                 library(stringr)
                                 library(readr)
                                 library(dplyr)
                                 library(tidylog)
                                 library(optparse)
                               })


option_list <- list(
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

opt <- parse_args(OptionParser(option_list = option_list))

files <- Sys.glob(paste0(opt$input, '/*.tsv'))

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
De_Novo_Junctions
chr
strand
Junctions_coords
Exons_coords
IR_coords
UCSC_LSV_Link"

read_majiq_out <- function(x) {
  read_tsv(
    x,
    col_names = strsplit(col_names, '\n')[[1]],
    comment = '#',
    cols(
      .default = col_character(),
      A5SS = col_logical(),
      A3SS = col_logical(),
      ES = col_logical(),
      Num_Junctions = col_double(),
      Num_Exons = col_double()
    )
  )
}

res <- lapply(files, read_majiq_out)

names(res) <- gsub(
  x = files,
  replacement = '\\1',
  pattern = file.path(opt$input, '(.*)_voila.tsv')
)

res <- bind_rows(res, .id = 'comparison') %>%
  separate_rows(
    "Junctions_coords",
    "E_dPSI_per_LSV_junction",
    "P_dPSI_beq_per_LSV_junction",
    "P_dPSI_leq_per_LSV_junction",
    "ref_E_PSI",
    "alt_E_PSI",
    sep = ';',
    convert = T
  )
# rank SJ by usage proportion
ref_rank <- res %>%
  group_by(comparison, LSV_ID) %>%
  group_modify(~dense_rank(.$ref_E_PSI) %>%
    tibble::enframe(value = 'ref_rank'))

alt_rank <- res %>%
  group_by(comparison, LSV_ID) %>%
  group_modify(~dense_rank(.$alt_E_PSI) %>%
    tibble::enframe(value = 'alt_rank'))

ranks <- alt_rank %>% inner_join(ref_ranke)

junction_pattern <- "(\\d+)-(\\d+)"
junctions_coords <- str_match(
  res$Junctions_coords, junction_pattern)[, c(2, 3)]

res['start'] <- junctions_coords[, 1]
res['end'] <- junctions_coords[, 2]
res <- res %>% select(
  comparison,
  chr,
  start,
  end.
  strand,
  P_dPSI_beq_per_LSV_junction,
  P_dPSI_leq_per_LSV_junction,
  E_dPSI_per_LSV_junction,
  LSV_ID)

write_csv(res, opt$output)
