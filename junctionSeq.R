# Title     : TODO
# Objective : TODO
# Created by: tbrittoborges
# Created on: 25/04/2018
suppressPackageStartupMessages(library(JunctionSeq));
options(stringsAsFactors = FALSE);

args <- commandArgs(TRUE)
threads <- as.integer(args[1])

decoder <- read.table(
  'junctionseq/decoder.tab', header=T, stringsAsFactors=F);
sample.files <-  paste0(
  "junctionseq/mergedOutput/",
  decoder$sample.ID,
  "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
message("Starting runJunctionSeqAnalyses (",date(),")");
jscs <- runJunctionSeqAnalyses(
    sample.files = sample.files ,
    sample.names = decoder$sample.ID,
    condition = decoder$group.ID,
    nCores = threads,
    verbose=TRUE,
    debug.mode = TRUE,
    flat.gff.file='junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz',
    # flat.gff.file='test/withNovel.forJunctionSeq.gff.gz',
    use.multigene.aggregates = TRUE,
    analysis.type = "junctionsAndExons");

message("Done with analysis (",date(),")");

dir.create(file.path('junctionseq/analysis/'))
message("Starting writeCompleteResults (",date(),")");
writeCompleteResults(
  jscs,
  outfile.prefix="junctionseq/analysis/",
  save.jscs = TRUE);
message("Done with writeCompleteResults (",date(),")");

dir.create(file.path('junctionseq/results/'))
message("Starting buildAllPlots (",date(),")");
bap <- buildAllPlots(
            jscs=jscs,
            FDR.threshold = 0.01,
            outfile.prefix = "junctionseq/results/",
            variance.plot = TRUE,
            ma.plot = TRUE,
            rawCounts.plot=TRUE,
            verbose = TRUE);
message("Done with  buildAllPlots (",date(),")");
