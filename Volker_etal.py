#!/usr/bin/env python3
# -*- coding: utf-8
"""
Created on 15:30 30/01/2018 2018 
Pipeline to process Majiq output table
.. usage:
python Volker_etal.py --file /Volumes/prj/Niels_Gehring/revision/RNPS1_Luc_workflow/RNPS1_Luc.deltapsi_deltapsi.tsv  --names  /Volumes/prj/Niels_Gehring/revision/RNPS1_Luc_workflow/luc.combined.luc.gtf.refmap

install maxentpy with pip install git+git://github.com/kepbod/maxentpy
"""
import os

import click as click
import numpy as np
import pandas as pd
import sys
from maxentpy import maxent_fast
import pybedtools


def load_gff_to_df(filename):
    header = ('NAME', 'SOURCE', 'TYPE', 'START', 'END',
              'SCORE', 'STRAND', 'FRAME', 'GROUP', 'empty')
    dtypes = {'END': int, 'FRAME': str, 'GROUP': str, 'ID': str, 'NAME': str,
              'NOTE': str, 'SCORE': str, 'SOURCE': str, 'START': int,
              'STRAND': str, 'TYPE': str, 'empty': str, }
    return pd.read_table(filename, names=header, na_values='.', dtype=dtypes,
                         comment='#')


def gff_query_junction(df, jun, chr=None, type='exon', strand='+'):
    q = 'NAME == "{chr}" & STRAND == "{strand}" & TYPE == "{type}" & START < ' \
        '{jun} & END > {jun} '
    return df.query(
        q.format(chr=chr, strand=strand, jun=jun, type=type))


def a5ss_exon_from_stringtie(row, dataset):
    if row['strand'] == '+':
        junction = row['junction_start']
    else:
        junction = row['junction_end']
    query = 'NAME == "{chr}" & STRAND == "{strand}" & TYPE == "exon" ' \
            '& START < {jun} &  {jun} < END '.format(
        chr=row['chr'], strand=row['strand'], jun=junction)
    result = dataset.query(query)
    if result.shape[0] != 1:
        return 1, 2
    return result.iloc[0, [3, 4]]


def explode(dataframe, columns, fill_value=''):
    """Decovolutes a the `columns` of a DataFrame

    .. note:
    https://stackoverflow.com/questions/45846765/efficient-way-to-unnest-
    explode-multiple-list-columns-in-a-pandas-dataframe
    :param dataframe:
    :param columns:
    :param fill_value:
    :return pandas.DataFrame:
    """
    if columns and not isinstance(columns, list):
        columns = [columns]
    # all columns except `lst_cols`
    idx_cols = dataframe.columns.difference(columns)

    # calculate lengths of lists
    lens = dataframe[columns[0]].str.len()

    if (lens > 0).all():
        # ALL lists in cells aren't empty
        return pd.DataFrame({
            col: np.repeat(dataframe[col].values,
                           dataframe[columns[0]].str.len())
            for col in idx_cols
        }).assign(**{col: np.concatenate(
            dataframe[col].values) for col in columns}).loc[:,
               dataframe.columns]
    else:
        # at least one list in cells is empty
        return pd.DataFrame({
            col: np.repeat(dataframe[col].values,
                           dataframe[columns[0]].str.len())
            for col in idx_cols
        }).assign(
            **{col: np.concatenate(dataframe[col].values) for col in columns}) \
                   .append(dataframe.loc[lens == 0, idx_cols]).fillna(
            fill_value) \
                   .loc[:, dataframe.columns]


def junction_type(row):
    """Junction type classification"""
    type_ = []

    j_start, j_end, exons, strand = row.loc[[
        'junction_start', 'junction_end', 'new_exons', 'strand']]
    if j_end - j_start <= 1:
        return ['IR']
    if not exons:
        return ['EMPTY']

    for i, (e_start, e_end) in enumerate(
            exon.split('-') for exon in exons):
        e_start, e_end = int(e_start), int(e_end)
        if j_start in (e_start, e_end):
            type_.append('JS')
            continue

        elif j_end in (e_start, e_end):
            type_.append('JE')
            continue

        if strand == '+':
            a = e_start - j_start
            b = j_end - e_end
            c = e_end - j_start
            d = j_end - e_start
        else:
            a = j_end - e_end
            b = e_start - j_start
            c = j_end - e_start
            d = e_end - j_start

        if a > 0 and b > 0 and c > 0 and d > 0:
            type_.append("ES")
        elif a < 0 and b > 0 and c > 0 and d > 0:
            type_.append("A5SS")
        elif a > 0 and b < 0 and c > 0 and d > 0:
            type_.append("A3SS")
        elif a < 0 and b < 0 and c > 0 and d > 0:
            type_.append("EX")
        else:
            type_.append("NONE")

    return type_


def check_exon_inclusion(group):
    canonical_idx = group['E(dPSI) per LSV junction'].idxmax()
    putative_EI_idx = group['E(dPSI) per LSV junction'].idxmin()

    if 'ES' not in group.loc[canonical_idx, 'junction_type']:
        # if there is no ES, nothing to look at
        return False
    elif 'IR' in group.loc[putative_EI_idx, 'junction_type']:
        # also, if there is IR in the canonical
        return False

    c_start, c_end = group.loc[
        canonical_idx, ['junction_start', 'junction_end']]
    j_start, j_end = group.loc[
        putative_EI_idx, ['junction_start', 'junction_end']]

    if c_end - c_start > j_end - j_start:
        return True
    else:
        return False


def bed_from_df(df, columns):
    # look at gtf to fasta from tophat
    # colum order gene, chr, strand, start, end

    bed = ''
    for idx, row in df.loc[:, columns].iterrows():
        name, chr, strand, start, end = row
        if strand == '+':
            bed = '\n'.join(
                [bed, f'{chr}\t{start - 3}\t{start + 6}\t{name}\t1\t{strand}'])

        elif strand == '-':
            bed = '\n'.join(
                [bed, f'{chr}\t{end - 7}\t{end + 2}\t{name}\t1\t{strand}'])
    return bed


def get_sequences(bed):
    pybed_obj = pybedtools.BedTool(bed, from_string=True)
    pybed_obj.sequence(
        fi='/Users/tbrittoborges/GRCh38_90.fa',
        name=True, s=True)

    with open(pybed_obj.seqfn) as f:
        seq = f.readlines()

    return [x.rstrip() for x in seq if not x.startswith('>')]


def split_exon(exon):
    if isinstance(exon, str):
        return (int(x) for x in exon.split('-'))
    return np.nan, np.nan


def first_event_exon_coord(row, event='A5SS', exon_label='new_exons'):
    if event not in row['junction_type']:
        return '0-0'  # None

    idx = row['junction_type'].index(event)
    return row[exon_label][idx]


def event_distance(row, event='A5SS'):
    try:
        exon_s, exon_e = split_exon(first_event_exon_coord(row, event))
    except AttributeError:  # not string
        return np.nan

    if row['strand'] == '+':
        return row['junction_start'] - exon_s
    else:
        return exon_e - row['junction_end']


@click.command()
@click.option('--file', type=click.Path(exists=True),
              help='Majiq output (tsv file)')
@click.option('--posterior_prob_filter', type=float, default=0.90,
              help='Filter junctions with P(dPSI) < {default}')
@click.option('--fold_change_filter', type=float, default=0.10,
              help='Filter junctions with |E(dPSI)| < {default}')
@click.option('--refmap_path', type=click.Path(exists=True),
              help='gtf.refmap file from StringTie, which contain the'
                   'mapping between transfrags and the reference'
                   'gene name, gene id and also the mapping class code')
@click.option('--gff_path', type=click.Path(exists=True),
              help='Path to Majiq gff file used for exon definition')
def main(file, posterior_prob_filter, fold_change_filter, refmap_path,
         gff_path):
    base_name = os.path.basename(file).split('.')[0]
    table = pd.read_table(file)
    dPSI_col = table.columns[table.columns.str.startswith('P(|dPSI|')][0]
    PSI_cols = table.columns[table.columns.str.endswith('E(PSI)')]

    print('Majiq output {} loaded with {} rows and {} columns.'.format(
        file, *table.shape))
    if refmap_path:
        # Use StringTie to fetch the Gene names and class code for the mapping
        names_table = pd.read_table(refmap_path)
        names_table['qry_id_list'] = names_table['qry_id_list'].str.split('|')
        names_table = explode(names_table, ['qry_id_list'])
        names_table.rename(columns={'qry_id_list': '#Gene Name'}, inplace=True)
        names_table.set_index('#Gene Name', inplace=True)
        table = table.join(names_table, on='#Gene Name')

    # process the convoluted columns
    table['Exons coords'] = table['Exons coords'].str.split(';')
    cols = ['Junctions coords', 'E(dPSI) per LSV junction', dPSI_col] + list(
        PSI_cols)

    for col in cols:
        table[col] = table[col].str.split(';')

    # Deconvolute the LSV/row to junction/row
    table = explode(table, cols)

    for col in cols[1:]:
        table[col] = table[col].astype(float)

    table['junction_start'] = (table['Junctions coords']
        .str.split('-', expand=True)[0]
        .astype(int))

    table['junction_end'] = (table['Junctions coords']
        .str.split('-', expand=True)[1]
        .astype(int))

    # Filter unprobable and low fold change events
    table['filtered'] = False
    table.loc[table['E(dPSI) per LSV junction'] > -fold_change_filter,
              'filtered'] = True
    table.loc[table[dPSI_col] < posterior_prob_filter,
              'filtered'] = True
    print(table.filtered.value_counts())

    # extract exons from gff
    gff = (pybedtools.BedTool(fn=gff_path)
        .filter(lambda x: x[2] == 'exon')
        .sort()
        .saveas())
    # need to normalise the bedfile to 6 columns
    table['dot'] = '.'
    table['idx'] = range(table.shape[0])

    bed_cols = ['chr', 'junction_start', 'junction_end', 'idx', 'dot', 'strand']
    a_bed = pybedtools.BedTool.from_dataframe(table.loc[:, bed_cols])
    intersect = (a_bed.intersect(b=gff, wa=True, wb=True, s=True)
        .to_dataframe(names=range(15)))
    intersect['new_exons'] = intersect[9].astype(str) + "-" + intersect[
        10].astype(str)

    exons = intersect.groupby([3], sort=False)['new_exons'].apply(
        lambda x: list(set(x)))

    # exons['new_exons'] = exons['new_exons'].str.join(' ')
    table = table.merge(
        right=exons.to_frame(), left_on=['idx'], right_index=True)

    # Add junction type for each exon in the majiq output
    # A5SS, A3SS, ES, EI, IR, EX

    table['junction_type'] = table.apply(junction_type, axis=1)

    exon_inclusion = table.groupby('LSV ID').apply(check_exon_inclusion)
    table['EI'] = table['LSV ID'].apply(exon_inclusion.get)
    for ass_type in 'A5SS A3SS ES IR EX'.split():
        table[ass_type] = table['junction_type'].apply(
            lambda x: ass_type in x)

    # Deep analysis
    # hack to deal with lack of typed Nan
    es_table = table.loc[table['ES']].copy()
    # a5ss_table = table.loc[table['A5SS']].copy()
    es_table['es_first'] = es_table.apply(
        first_event_exon_coord, event='ES', axis=1)
    es_table['ES_start'] = (es_table['es_first']
        .str.split('-', expand=True)[0]
        .astype(int)
        .replace(0, np.nan))
    es_table['ES_end'] = (es_table['es_first']
        .str.split('-', expand=True)[1]
        .astype(int)
        .replace(0, np.nan))

    a5ss_table = table.loc[table['A5SS']].copy()

    a5ss_table['A5SS_first'] = a5ss_table.apply(first_event_exon_coord, axis=1)
    a5ss_table['A5SS_start'] = (a5ss_table['A5SS_first']
        .str.split('-', expand=True)[0]
        .astype(int)
        .replace(0, np.nan))
    a5ss_table['A5SS_end'] = (a5ss_table['A5SS_first']
        .str.split('-', expand=True)[1]
        .astype(int)
        .replace(0, np.nan))

    for t, component in zip([table, es_table, a5ss_table],
                            ['junction', 'ES', 'A5SS']):
        start, end = f'{component}_start', f'{component}_end'
        bed_ = bed_from_df(t, ['#Gene Name', 'chr', 'strand', start, end])
        t[f'{component}_seq'] = get_sequences(bed_)
        t[f'{component}_score5'] = t[f'{component}_seq'].apply(
            maxent_fast.score5)
        t['dinucleotide'] = t[f'{component}_seq'].str[3:5]

    a5ss_table['A5SS_dist'] = table.apply(event_distance, axis=1)

    # output junctions counts per gene
    cols = {'junction': ['chr', 'junction_start', 'junction_end', 'strand',
                         'junction_type', 'Exons coords', 'EI', 'LSV ID',
                         'new_exons',
                         'junction_seq', 'dinucleotide', 'junction_score5',
                         '#Gene Name', 'Gene ID', 'filtered', 'A5SS', 'A3SS',
                         'ES', 'IR', 'EX', 'E(dPSI) per LSV junction', dPSI_col,
                         *PSI_cols],

            'ES': ['chr', 'ES_start', 'ES_end', 'strand', 'junction_start',
                   'junction_end', 'Exons coords', 'EI', 'LSV ID', 'ES_seq',
                   'ES_score5', 'filtered', 'E(dPSI) per LSV junction',
                   'new_exons', dPSI_col, *PSI_cols],

            'A5SS': ['chr', 'A5SS_start', 'A5SS_end', 'strand', 'new_exons',
                     'junction_start', 'junction_end', 'Exons coords',
                     'EI', 'LSV ID', 'A5SS_seq', 'A5SS_score5',
                     'A5SS_dist', 'filtered', 'E(dPSI) per LSV junction',
                     dPSI_col, *PSI_cols]}

    if refmap_path:
        table.loc[table.filtered, 'ref_gene_id'].value_counts().to_csv(
            '{}.junction_per_genes.csv'.format(base_name))
        for col in cols:
            cols[col].extend('ref_gene_id ref_id class_code'.split())

    for t, component in zip([table, es_table, ],  # a5ss_table
                            ['junction', 'ES', ]):  # 'A5SS'
        t = t.loc[:, cols[component]].copy()
        t.columns = t.columns.str.replace('per LSV junction', '')
        t.to_csv('{}_{}.csv'.format(base_name, component),
                 index=False)


if __name__ == "__main__":
    main()
