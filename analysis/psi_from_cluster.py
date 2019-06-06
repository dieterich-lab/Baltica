#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 14:20 2019-06-06 2019 

"""

import pandas as pd

lc = pd.read_csv()
ju = pd.read_csv()

lc = lc[
    (lc['p.adjust'] < 0.05) &
    (lc['deltapsi'].abs() > 0.1)]

lc.end = lc.end - 1
lc.genes[lc.genes.isna()] = ''

results = {}

for group, x in lc.groupby(['contrast', 'cluster']):
    comp = group[0]
    psi = {}
    ctr, trt = comp.split('-vs-')

    if x.shape[1] == 1:
        continue

    ref_j = x.loc[x[ctr].idxmax()].copy()
    x = x.drop(x[ctr].idxmax())
    if x.empty:
        continue

    overlap = (x.start >= ref_j.start) & (x.end <= ref_j.end)
    x = x.drop(overlap.index[overlap == False])
    if x.empty:
        continue
    alt_j = x.loc[x[trt].idxmax()].copy()

    ref_counts = ju.query('chr == @ref_j.chr & start == @ref_j.start & end == @ref_j.end')

    alt_counts = ju.query('chr == @alt_j.chr & start == @alt_j.start & end == @alt_j.end')

    results[group] = {
        'gene': ref_j.genes,
        'ref_chr': ref_j.chr,
        'ref_start': ref_j.start,
        'ref_end': ref_j.end,
        'ref_ctr_counts': float(
            ref_counts.loc[ref_counts.index, ref_counts.columns.str.startswith(ctr)].mean(1)),
        'alt_ctr_counts': float(
            alt_counts.loc[alt_counts.index, alt_counts.columns.str.startswith(ctr)].mean(1)),
        'alt_chr': alt_j.chr,
        'alt_start': alt_j.start,
        'alt_end': alt_j.end,
        'ref_trt_counts': float(
            ref_counts.loc[ref_counts.index, ref_counts.columns.str.startswith(trt)].mean(1)),
        'alt_trt_counts': float(
            alt_counts.loc[alt_counts.index, alt_counts.columns.str.startswith(trt)].mean(1)),
    }

results = pd.DataFrame.from_dict(results).T

results['ctr_psi'] = results['alt_ctr_counts'] / (
            results['ref_ctr_counts'] + results['alt_ctr_counts'])
results['trt_psi'] = results['alt_trt_counts'] / (
            results['ref_trt_counts'] + results['alt_trt_counts'])

