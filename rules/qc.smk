# -*- coding: utf-8
"""
Created on 12:00 07/01/2020
Snakemake file for quality control of RNA-Sequencing dataset focused on junction quality.
"""

__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Thiago Britto Borges"
__email__ = "tbrittoborges@uni-heidelberg.de"
__license__ = "MIT"

rule all:
    input:
        expand("qc/fastqc/{sample}_fastqc.zip", sample=sample),
        "ref.bed12",
        "qc/rseqc/gene_body_coverage.pdf",
        "qc/rseqc/gene_body_coverage.pdf",
        "qc/rseqc/inner_distance_{sample}.pdf",
        "qc/rseqc/{sample}.GC_plot.pdf",
        "qc/rseqc/{sample}.DupRate_plot.pdf",
        "qc/rseqc/{sample}.splice_junction.pdf",
        "qc/rseqc/{sample}.splice_events.pdf",
        "qc/rseqc/{sample}.junctionSaturation_plot.pdf",
        "qc/rseqc/{sample}_infer_experiment.txt",
        "qc/rseqc/{sample}.bam_stat.txt",
        "qc/resec/{sample}.read_distribution.txt"

# from https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/markduplicates.html
# rule mark_duplicates:
#     input:
#         "mapped/{sample}.bam"
#     output:
#         bam="dedup/{sample}.bam",
#         metrics="dedup/{sample}.metrics.txt"
#     log:
#         "logs/picard/dedup/{sample}.log"
#     params:
#         "REMOVE_DUPLICATES=true"
#     wrapper:
#         "0.45.1/bio/picard/markduplicates"

# java -Xmx8G -jar /biosw/picard-tools/2.5.0/picard.jar MarkDuplicates INPUT={input} OUTPUT={output.bam}  M={output.metrics}


# from https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastqc.html
rule fastqc:
    input:
        "reads/{sample}.fastq"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.45.1/bio/fastqc"

# conda install ucsc-gtftogenepred
# conda install ucsc-genepredtobed
rule ref_annotation_gtf_to_bed:
    input:
        config.ref
    output:
        "ref.bed12"
    shell:
        """
        wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
        wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
        chmod 755 gtfToGenePred genePredToBed
        ./gtfToGenePred {input} temp.genepred
        ./genePredToBed temp.genepred {output}
        rm temp.genepred gtfToGenePred genePredToBed
        """

rule rseqc_gene_body_coverage:
    input:
        bed = "ref.bed12",
        bam = ','.join(expand("mapped/{sample}.bam", sample=sample))
    output:
        "qc/rseqc/gene_body_coverage.pdf"
    shell:
        """
        geneBody_coverage.py -r {input.bed} -i {input.bam} -o {output}
        """

rule rseqc_inner_distance:
    input:
        ref = "ref.bed12",
        bam = "mapped/{sample}.bam"
    output: "qc/rseqc/{sample}.inner_distance_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        inner_distance.py -i {input.bam} -o {params.prefix} -r {input.ref}
        """

rule rseqc_read_gc:
    input: "mapped/{sample}.bam"
    output: "qc/rseqc/{sample}.GC_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        read_GC.py -i {input} -o {params.prefix}
        """

### fix above output
rule rseqc_read_duplication:
    input: "mapped/{sample}.bam"
    output: "qc/rseqc/{sample}.DupRate_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        read_duplication.py -i {input} -o {params.prefix}
        """


rule rseqc_junction_annotation:
    input:
        ref = "ref.bed12",
        bam = "mapped/{sample}.bam"
    output:
        "qc/rseqc/{sample}.splice_junction.pdf",
        "qc/rseqc/{sample}.splice_events.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        junction_annotation.py -i {input.bam} -r {input.bed} -o {output}
        """

rule rseqc_junction_saturation:
    input:
        ref = "ref.bed12",
        bam = "mapped/{sample}.bam"
    output: "qc/rseqc/{sample}.junctionSaturation_plot.pdf"
    params:
        prefix = "qc/rseqc/{sample}"
    shell:
        """
        junction_saturation.py -i {input.bam} -r {input.bed} -o {params.prefix}
        """

rule rseqc_infer_experiment:
    input:
        ref = "ref.bed12",
        bam = "mapped/{sample}.bam"
    output: "qc/rseqc/{sample}_infer_experiment.txt"
    shell:
        """
        infer_experiment.py -r {input.ref} -i {input.bam} > {output}
        """

rule rseqc_bam_stat:
    input: "mapped/{sample}.bam"
    output: "qc/rseqc/{sample}.bam_stat.txt"
    shell:
        """
        bam_stat.py -i {input} > {output}
        """


rule read_distribution:
    input:
        bam = "mappings/{sample}.bam",
        ref = "ref.bed12"
    output:
        "qc/rseqc/{sample}.read_distribution.txt"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.bed} &> {output}
        """

# rule multiqc:
#     input:

#     output:
#     shell:
#         """
#         multiqc --fullnames -d --dirs-depth 4 fastqc_new_fastq/ fastqc_mappings/ /prj/Niels_Gehring/newData_March_2018/workflow_spike_ins/mapping/ /prj/Niels_Gehring/newData_July_2019_CASC3/workflow/mapping/  /prj/Niels_Gehring/newData_July_2019_CASC3/workflow/rrna_free_reads/  /prj/Niels_Gehring/newData_March_2018/workflow_spike_ins/rrna_free_reads/ --force --verbose  --ignore "MATE"
#         """
