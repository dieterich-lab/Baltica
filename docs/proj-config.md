# Baltica project configuration:


Baltica requires a project configuration file as input. For a template see [here](https://raw.githubusercontent.com/dieterich-lab/Baltica/master/baltica/config.yml). For a programmatic solution to generate the configuration, use the script [baltica/write_new_config.R](https://github.com/dieterich-lab/Baltica/blob/8bc5fe5f71e948b3971ea5db5f1456b2d9e2f838/baltica/write_new_config.R). Method specific parameters are detailed in the [Workflow implementation](workflows.md) page.

!!! warning
   Baltica project requires that configuration been update. Users should use the full path to files. The sample names for the `samples` parameters uses a underscore (\_) to separate the sample name and replicate number, like in mix\_1, so the sample should not contain spaces or underscores. 

!!! note
    The `Required` column flags parameters without default.

## General parameters description

| Parameter | Description | Required |
|---|---|---|
| path | project path | &#10003; |
| sample_path | path to the parent directory for aligment files | &#10003; |
| samples | sample name and directory formated as "{sample_name}:   {path_to_sample) | &#10003; |
| contrasts | list of contrasts names ([format](proj-config.md#pairwise-comparisons)) | &#10003; |
| assembly | assembly name on UCSC Browser |  |
| strandness | one of fr-firststrand, fr-secondstrand or none (unstranded) | &#10003; |
| read_len | maximun read length  | &#10003; |
| ref | path to reference annotation in the GTF format | &#10003; |
| ref_fa | path to reference annotation in the FASTA format | &#10003; |
| project_authors | project author name, used in the report |  |
| project_title | project title name, used in file names and report |  |
| orthogonal_result | result from Nanopore-seq in GFF or BED with a valid score column, and optionally a comparisons column with   contrasts |  |


## Pairwise comparisons

```yaml
contrasts:
  {case1}-vs-{control}:
    - {case1}
    - {control} 
  {cas2}-vs-{control}:
    - {case2}
    - {control}
```

!!! note
    junctionseq and leafcutter support more complex experimental designs, which were not implemented in Baltica.


## majiq specific paramers

[majiq manual](https://biociphers.bitbucket.io/majiq/MAJIQ.html#builder)

| Parameter | Original parameter | Description |
|---|---|---|
| __rule majiq_build:__ (`majiq deltapsi`)| | |
| majiq_min_experiments | --min-experiments | minimum number of experiments to filter with --minreads (default 1.0 - all experiments) |
| majiq_minreads | --minreads | Discard SJ with less than `--minreads` reads |
| majiq_min_denovo | --min-denovo | Discard novel SJ with less than `--min-denovo` reads |
| __rule majiq_deltapsi__ (`majiq deltapsi`) | | |
| majiq_minreads | --minreads | same as above |
| __rule majiq_voila__ (`voila tsv`) | | |
| majiq_non_changing_threshold | --non-changing-threshold | |
| majiq_threshold | --threshold ||


## junctionseq specific paramers
[qorts manual](https://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf)
[junctionseq manual](http://hartleys.github.io/JunctionSeq/doc/JunctionSeq.pdf)

| Parameter | Original parameter | Description |
|---|---|---|
| __rule junctionseq_qc__ (`qorts QC`) | | |
| is_single_end | --singleEnded | Flag single end libraries |
| __rule junctionseq_merge__ (`qorts mergeNovelSplices`) | | |
| junctionseq_mincount | --minCount | Discard SJ with less than `--minCount` reads |

## leafcutter specific paramers

[regtools manual](https://regtools.readthedocs.io/en/latest/)
[leafcutter manual](http://davidaknowles.github.io/leafcutter/articles/Usage.html)

| Parameter | Original parameter | Description |
|---|---|---|
| __rule leafcutter_bam2junc__ (`regtools junctions extract`) | | |
| leafcutter_minimum_anchor_length | -a | Discard reads with overanging length lower than `-a` |
| leafcutter_minimum_intron_size | -i | Minimum intron size |
| leafcutter_maximum_intron_size | -I | Maximum intron size |
| __rule leafcutter_differential_splicing__ | | |
| leafcutter_min_coverage | --min_coverage | desc |
| leafcutter_min_samples_per_group | -g | Discard SJ used in less than `-g` samples | |
| leafcutter_min_samples_per_intron | -i | Discard SJ with less than `--min_coverage` in `-i` samples in each group |

## rmats specific paramers

[rmats manual](https://github.com/Xinglab/rmats-turbo#all-arguments)

| Parameter | original parameter | Description |
|---|---|---|
| __rule rmats_run__ (`rmats.py`) | | |
| rmats_allow_clipping | --allow-clipping | Allow clipped reads |
| rmats_variable_read_length | --variable-read-length | Allow reads with variable-length  |
| rmats_novel_ss | --novelSS | Allow the detecting unannotated SJ |
| rmats_extra | none | Pass extra arguments to `rmats.py` |
