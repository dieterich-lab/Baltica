source('../baltica/utils.R')

get_changing_introns <- function() {
  gtf=rtracklayer::import.gff2('/beegfs/homes/tbrittoborges/Baltica/data/SIRV.gtf')
  # > length(unique(gtf$transcript_id))
  # [1] 101
  # > length(unique(gtf$gene_id))
  # [1] 7
  changing=read.table('/prj/Niels_Gehring/baltica_benchmark/Christoph/Changing_SIRV.txt')
  # nrow(changing)
  # 52
  # if a transcript changes it changes, it is different in the three mixes
  # some transcripts do not change
  # some are absent in the mixes, but present in the annotation
  # all changes in molarity should be detected, the smallest change is from 1M to 0.5M
  ex_tx=filter_multi_exon(gtf)
  introns=get_introns(ex_tx)
  # > length(introns)
  # [1] 482
  # > length(unique(introns))
  # [1] 138
  mcols(introns)$changing=names(introns) %in% changing$V1
  hits=findOverlaps(introns, type='equal', select='first')
  unique_names=split(names(introns), hits)
  unique_changing=split(introns$changing, hits)
  introns=introns[unique(hits)]
  mcols(introns)$changing=unlist(lapply(unique_changing, any))
  # > table(introns$changing)
  # FALSE  TRUE
  #    40    98
  introns=unname(introns)
  expand=expand.grid(seq_along(introns), c("mix2-vs-mix1", "mix3-vs-mix1", "mix3-vs-mix2"))
  introns=introns[expand$Var1]
  mcols(introns)$comparison = expand$Var2
  introns
  
}

introns <- get_changing_introns()
mcols(introns)$score <- ifelse(mcols(introns)$changing, 1, 0)
rtracklayer::export.gff2(introns, '~/Baltica/data/SIRV_changing_introns.gtf')
