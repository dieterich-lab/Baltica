#!/usr/bin/env Rscript --vanila
# Title : JunctionSeq analysis
# Objective : Call functions for JunctionSeq analysis in context of Baltica
# Created by: tbrittoborges
# Created on: 25/04/2018
suppressPackageStartupMessages(library(JunctionSeq))
options(stringsAsFactors = FALSE);

args <- commandArgs(TRUE)
threads <- snakemake@config[["threads"]]

decoder <- read.table(
  'junctionseq/decoder.tab', header = T, stringsAsFactors = F);

config <- yaml::read_yaml(snakemake@params[['configpath']])

sample.files <- paste0(
  "junctionseq/mergedOutput/",
  decoder$sample.ID,
  "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
message("Starting runJunctionSeqAnalyses (", date(), ")");

for (x in names(config$contrasts)){
    message("Starting runJunctionSeqAnalyses (", date(), ") for ", x)

    subset_ <- decoder$group.ID %in% config$contrasts[[x]]

    jscs <- runJunctionSeqAnalyses(
      sample.files = sample.files[subset_],
      sample.names = decoder[subset_, 'sample.ID'],
      condition = decoder[subset_, 'group.ID'],
      nCores = threads,
      verbose = TRUE,
      debug.mode = TRUE,
      flat.gff.file = 'junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz',
      use.multigene.aggregates = TRUE,
      analysis.type = "junctionsOnly");

    out_path <- file.path('junctionseq/analysis/', x)
    dir.create(out_path)
    message("Starting writeCompleteResults (", date(), ") for ", x)
    writeCompleteResults(
      jscs,
      outfile.prefix = out_path,
      save.bedTracks = FALSE,
      save.jscs = TRUE)

    message("Done with writeCompleteResults for ", x, (", date(),  ")")

}
message("Done with all comparisons (", date(), ")")