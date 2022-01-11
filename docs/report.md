# Baltica output

Baltica framework produces two files as output:  
    - an R markdown report
    - an excel spreadsheet  

!!! note
    If available, the orthogonal dataset is treated as a new method named `orthogonal.`

## Baltica table spreadsheet

- `results/baltica_table_{proj_name}.xlsx`

The spreadsheet contains the complete set of coordinate output by methods and comparisons. In addition, there are a column for the combination of methods and comparisons plus the columns for the annotation:

- coordinates: junction genomic coordinate in the format: `{chr}:{start}-{end}:{strand}` (strand omitted if none)
- score columns: in the format: `{method}_{comparisons}`
- is_novel: whether the splicing junction is or not annotated
- gene_name: the gene name obtained from the de novo annotation workflow
- transcript_name: transcript name from the de novo annotation
- class_code: transcript class association to the reference annotation transcript, please see [Fig 1 in the GFF Utilities paper](https://f1000research.com/articles/9-304/v2) for details
- exon_number: pairs of exon numbers from the de novo annotation. First of the pair is the donor exon if the feature is the positive strand; otherwise acceptor
- as_type: type of AS for each junction exon skipping (ES), alternative 3' splice site (A3SS), alternative 5' splice site (A5SS)

Currently, the HTML report comprises two sections:

## Common splice junctions

The [upset plot](https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html#upset-plot) shows the combination of distinct sets of calls (score > 0.95) from each method and contrast. The plot helps to compare the common calls among sets. The complement sets are ignored, as these sets usually have a high number.

## Baltica table

This interactive HTML table provides the top 1,000 (or `baltica_max_table` in the configuration file) sorted by the sum of the scores. Extra annotation is available upon clicking on ▶. In addition, the coordinates columns link to the UCSC genome browser. Regional URL for UCSC GB can be selected with `ucsc_url`, and assembly should be selected with the `assembly` configuration

## Baltica report configuration

Change the following options on your project configuration to change the report:  

- project_authors: name of the persons running the project
- project_title: name of the files and report title
- baltica_max_table: maximum number of rows on the HTML table
- assembly: assembly used for linking with the genome browser
- ucsc_url: URL for the genome browser, like `http://genome-euro.ucsc.edu` for the European mirror

## Reproducibility

This section provides the information necessary to reproduce the report, including project configuration and R package version.  
