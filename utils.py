#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 14:59 19/01/2018 2018 

"""
from collections import namedtuple
from itertools import combinations

import numpy as np
import pandas as pd

GenomicRegion = namedtuple('GenomicRegion', ['chr', 'strand', 'start', 'end'])


def create_mapping(config):
    """
    Generate the mapping bettwen samples and replicates 
    :param config: snakemake configuration file
    :return dict: mapping bettwen sample and replicates
    """
    names = config["samples"].keys()
    conditions = sorted(set([x.split('_')[0] for x in names]))

    return {c: [x for x in names if x.split('_')[0] == c] for c in conditions}


def all_against_all(config):
    """

    :param config: 
    :return: 
    """
    names = config["samples"].keys()
    conditions = sorted(set([x.split('_')[0] for x in names]))
    return ['{}_{}'.format(*x) for x in combinations(conditions, 2)]


def region(chr, strand='+', start=None, end=None):
    """
    Validate the datatype for each value to facilitate DataFrame query mapping
    and facilitates the creation of open ended GenomicRegion objects

    :param str chr: chromossome
    :param str strand: plus or minus strand
    :param int or None start: position start or None
    :param int or None end:
    :return GenomicRegion:
    """
    return GenomicRegion(str(chr), strand, int(start), int(end))


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


def check_for_contigous_junctions(
        data, chr, strand, start, end):
    """Assert that an exon has support of reads that overlap introns,
    aka, reads from the BallGown i_data file"""

    start_support = not data.query(
        f'chr == "{chr}" & strand == "{strand}" & end == {start - 1}').empty
    end_support = not data.query(
        f'chr == "{chr}" & strand == "{strand}" & start == {end + 1}').empty

    return start_support and end_support


def exon_within(data, chr, start, end):
    """Query for Exons within a given genomic coordinate"""
    query = f'chr == "{chr}" & start > {start} & end < {end}'

    result = data.query(query)
    # simplest case:
    # if the exon overlaps with other multiples exons,
    # select the one with highest mrcounts
    if not result.empty:
        #         return result.loc[result['mrcount'].idxmax(),
        #                           ['chr', 'start', 'end']].tolist()
        return result


def exon_contains(data, chr, position, boundary=False):
    """Query Exons that contains a given genomic coordinate"""
    if boundary:
        query = f'chr == "{chr}" & start <= {position} <= end'
    else:
        query = f'chr == "{chr}" & start < {position} < end'

    return data.query(query)


def all_exons(data, chr, start, end, introns_data=None):
    """Query all elements involved in a junction"""
    df = pd.concat([
        exon_contains(data, chr, start, boundary=True),
        exon_within(data, chr, start, end),
        exon_contains(data, chr, end, boundary=True)
    ])

    df['junction_support'] = [check_for_contigous_junctions(
        introns_data, row['chr'], row['strand'], row['start'], row['end'])
        for _, row in df.iterrows()]

    if df.empty:
        return df

    #  ignore support for first and last exons
    # this is an approximation, this should be done
    # check 5p junction support for the first exon
    # check 3p junction support for the last exon

    df.loc[df['start'].idxmin(), 'junction_support'] = True
    df.loc[df['end'].idxmax(), 'junction_support'] = True

    return df
