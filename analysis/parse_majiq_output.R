#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "majiq/voila/*_voila.tsv",
    help = "Path with glob character to Majiq result files. [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "majiq/majiq_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    type = "double",
    default = 0.90,
    help = "Discard junction with probability threshold < than --cutoff [default %default]"
  )
)
# enabling both cli or snakemake input
tryCatch(
  {
    opt <- list(
      input = snakemake@input,
      output = snakemake@output[[1]],
      cutoff = snakemake@params[['cutoff']]
      )
    files <- opt$input
    file_names <- gsub(
     x = opt$input,
     pattern = 'majiq/voila/(.+)_voila.tsv',
     replacement = '\\1')

  }, error = function(e) {

    opt <- parse_args(OptionParser(option_list = option_list))
    files <- Sys.glob(opt$input)
    file_names <- gsub(
      x = files,
      replacement = '\\1',
      pattern = sub(x=opt$input, pattern='\\*', replacement = '(.*)')
)
})


# rename column names from majiq result due to the presence of spaces
read_majiq_out <- function(x) {
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

  read_tsv(
    x,
    col_names = strsplit(col_names, '\n')[[1]],
    comment = '#',
    skip = 1,
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

names(res) <- file_names

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
# flag canonical SJ
res <- res %>%
  arrange(comparison, LSV_ID, ref_E_PSI) %>%
  group_by(comparison, LSV_ID) %>%
  mutate(is_canonical = row_number() == 1) %>%
  ungroup()


junction_pattern <- "(\\d+)-(\\d+)"
junctions_coords <- str_match(
  res$Junctions_coords, junction_pattern)[, c(2, 3)]

res['start'] <- junctions_coords[, 1]
res['end'] <- junctions_coords[, 2]

message('Number of junctions output by Majiq ', nrow(res))
res <- dplyr::select(res, c(
  chr,
  start,
  end,
  comparison,
  P_dPSI_beq_per_LSV_junction,
  strand,
  LSV_ID,
  P_dPSI_leq_per_LSV_junction,
  E_dPSI_per_LSV_junction,
  ref_E_PSI,
  alt_E_PSI,
  is_canonical)
) %>%
  filter(P_dPSI_beq_per_LSV_junction > opt$cutoff) %>%
  mutate(method = 'majiq')

message('Number of junctions after filtering ', nrow(res))

write_csv(res, opt$output)