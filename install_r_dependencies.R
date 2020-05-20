#!/usr/bin/env Rscript

# Objective : Install R dependencies for Baltica analysis
# Created by: tbrittoborges
# Created on: 19.05.20
# Based on: https://github.com/dieterich-lab/circtools/blob/master/scripts/install_R_dependencies.R
args <- commandArgs(trailingOnly = TRUE)
cran <- "https://cran.uni-muenster.de/"
if (!is.na(args[1])){
  cran <- args[1]
}

options(
  repos = c(CRAN = cran))
Sys.setenv(R_INSTALL_STAGED = FALSE) # requires 3.6


install_many <- function (pkgs, install_fun, ...){
  pkgs <- pkgs[!pkgs %in% installed.packages()[,1L]]
  if (length(pkgs) > 0L)
    install_fun(pkgs, quietly = TRUE, ...)

}
pkgs.cran <- c(
  "UpSetR", "dplyr", "openxlsx", "optparse",  "stringr",  "tidyr", "readr"
)
# First install packages from CRAN
install_many(pkgs.cran, install.packages)

# Next, packages from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")
pkgs.bioc <- c(
  "GenomicRanges", "Rsamtools", "UpSetR", "yaml"
)
install_many(pkgs.bioc, BiocManager::install)
