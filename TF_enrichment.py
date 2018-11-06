#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 11:23 06/11/2018 2018 
Determining which motif(s) are relatively enriched in one set of promoters
compared to another set of promoters.
"""
import pathlib

import gffutils
import pandas as pd
import pybedtools
from gffutils.pybedtools_integration import tsses

selected_genes_path = '/Volumes/prj/SFB_OMICS_Mouse_disease_models/' \
                      'RNA/tbb_DE/DE_genes_list/'
selected_genes_patern = '*_selected_genes.csv'
database_path = '/Volumes/tbrittoborges/GRCm38.90.gtf.sql'
upstream_distance = 10000
downstream_distance = 100
chr_size_path = '/Volumes/tbrittoborges/GRCm38_90_chr_size.tab'
regulatory_feat_path = '/Volumes/prj/dorina2/raw/ensembl_regulatory/' \
                       'mus_musculus.GRCm38.heart_adult_8_weeks.' \
                       'Regulatory_Build.regulatory_activity.20180516.gff'
output_directory = '/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/tss/'
fout = '/Volumes/prj/SFB_OMICS_Mouse_disease_models/RNA/tbb_DE/tss/GRCm38' \
       '_90_{}_{}.fa'
fin = '/Volumes/biodb/genomes/mus_musculus/GRCm38_90/GRCm38_90.fa'

print('Loading transcript database')
db = gffutils.FeatureDB(
    database_path, keep_order=True)
tss = tsses(db)

tss = (tss.slop(
    l=upstream_distance,
    r=downstream_distance,
    g=chr_size_path)
       .sort()
       .saveas())

print('Loading regulatory dataset')
regulatory_feat = pybedtools.BedTool(
    regulatory_feat_path)
regulatory_feat = (regulatory_feat
                   .sort()
                   .filter(
    lambda tss: tss.attrs['activity'] == 'ACTIVE')
                   .saveas())

print('Loading regulatory dataset')
tss = tss.intersect(
    regulatory_feat,
    wa=True,
    sorted=True,
    stream=True).saveas()

print('Loading DE genes')
de_gene_p_condition = {}
path = pathlib.Path(selected_genes_path)
for file in path.glob(selected_genes_patern):
    cond = file.name.replace('_selected_genes.csv', '')
    de_gene_p_condition[cond] = set(pd.read_csv(file, index_col=0).index)

for cond in de_gene_p_condition:
    print(cond, ' number of DE genes', len(de_gene_p_condition[cond]))

# todo add max(padj) to gene name
print('processing promoter sequence')
for cond in de_gene_p_condition:
    de_seq = []
    control = []

    for row in tss:
        if row.attrs['gene_id'] in de_gene_p_condition[cond]:
            de_seq.append(row)
        else:
            control.append(row)

    pybedtools.BedTool(de_seq).sequence(fin).save_seqs(
        fout.format(cond, 'input'))
    pybedtools.BedTool(control).sequence(fin).save_seqs(
        fout.format(cond, 'control'))
