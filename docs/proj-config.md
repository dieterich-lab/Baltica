# Baltica project configuration:


Baltica requires a project configuration file as input. For a template see [here](https://raw.githubusercontent.com/dieterich-lab/Baltica/master/baltica/config.yml). For a programmatic solution to generate the configuration, use the script [baltica/write_new_config.R](https://github.com/dieterich-lab/Baltica/blob/8bc5fe5f71e948b3971ea5db5f1456b2d9e2f838/baltica/write_new_config.R). Method specific parameters are detailed in the [Workflow implementation](workflows.md) page.

!!! warning
    You are required to update configuration. Use full path to files.

!!! note
    The `Required` column flags parameters without default.

## Parameters description

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
| orthognal_result | result from Nanopore-seq in GFF or   BED with a valida score column, and optionally a comparisons column with   contrasts |  |


<!-- | **majiq specific parameters** |  |  |
| **junction specific parameters** |  |  |
| **stringtie specific parameters** |  |  |
| **leafcutter specific parameters** |  |  |
| leafcutter_min_samples_per_group |  |  |
| leafcutter_min_samples_per_intron |  |  |
| leafcutter_min_coverage |  |  |
| leafcutter_min_cluster_reads |  |  |
| leafcutter_max_intron_length |  |  |
| majiq_threshold |  |  |
| **majiq specific parameters** |  |  |
| majiq_minreads |  |  | -->

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
    junctionseq and leafcutter support more complex experimental designs, which were not yet implemented in Baltica.
