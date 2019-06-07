library(rtracklayer)

gtf <- rtracklayer::import( '/biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf' )
exons <- as.data.frame(gtf[gtf$type == 'exon', ])

write.csv(
  exons[, c(1, 2, 3, 5, 10, 12, 17, 26)],
  file= 'Thiago/Baltica/exons.csv', )


# We compress it deleting the .csv
system("gzip Thiago/Baltica/exons.csv")
