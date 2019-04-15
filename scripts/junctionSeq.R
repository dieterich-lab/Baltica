# Title     : TODO
# Objective : TODO
# Created by: tbrittoborges
# Created on: 25/04/2018
BiocManager::install('JunctionSeq', version = '3.8')

suppressPackageStartupMessages(library(JunctionSeq));
options(stringsAsFactors = FALSE);

threads <- as.integer(10)
setwd('/prj/Andre_Schneider_BNauheim/raw_silke_kreher/newBatchMar2019/Baltica')
decoder <- read.table(
  'junctionseq/decoder.tab', header=T, stringsAsFactors=F);
sample.files <-  paste0(
  "junctionseq/mergedOutput/",
  decoder$sample.ID,
  "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
message("Starting runJunctionSeqAnalyses (",date(),")");

if ('WT' %in% decoder$group.ID ) {
  decoder$group.ID <- relevel(
    as.factor(decoder$group.ID),
    'WT')
}

jscs <- runJunctionSeqAnalyses(
  sample.files = sample.files ,
  sample.names = decoder$sample.ID,
  condition = decoder$group.ID,
  nCores = threads,
  verbose=TRUE,
  debug.mode = TRUE,
  flat.gff.file='junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz',
  use.multigene.aggregates = TRUE,
  analysis.type = "junctionsOnly");

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
try(
  buildAllPlots(
    jscs=jscs,
    FDR.threshold = 0.05,
    outfile.prefix = "junctionseq/results/",
    variance.plot = TRUE,
    ma.plot = TRUE,
    rawCounts.plot=TRUE,
    verbose = TRUE)
)
message("Done with  buildAllPlots (",date(),")");