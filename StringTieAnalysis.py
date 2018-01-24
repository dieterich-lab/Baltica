#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 15:00 19/01/2018 2018 

e_data.ctab ->
i_data.ctab ->
"""
def check_for_contigous_junctions(
        data, chr, strand, start, end):
    """Assert that an exon has support of reads that overlap introns,
    aka, reads from the BallGown i_data file"""

    start_support = not data.query(
        f'chr == "{chr}" & strand == "{strand}" & end == {start - 1}').empty
    end_support = not data.query(
        f'chr == "{chr}" & strand == "{strand}" & start == {end + 1}').empty

    return start_support and end_support


def query_df(data, query):
    """

    :param pandas.Dataframe data:
    :param str query:
    :return:
    """
    return data.query(query)


def junction_support(intron_data, genomic_region):
    """
    Check junction support for a exon start or end.

    :param pandas.DataFrame intron_data: dataframe with the intro data
    :param genomic_region: the exon's genomic region
    :return bool: whether the region is supported by junctions
    """
    if genomic_region.start and not genomic_region.end :
        query = f'''chr == "{chr}" & strand == "{strand}" 
        & end == {start - 1}'''.format(genomic_region._asdict())

    elif not genomic_region.start and genomic_region.end:
        query = f'''chr == "{chr}" & strand == "{strand}"
         & end == {start - 1}'''.format(genomic_region._asdict())

    else:
        raise NotImplemented
        # return junction_support(intron_data, genomic_end) and # return
        # junction_support(intron_data, genomic_start)

    return not query_df(intron_data, query)
