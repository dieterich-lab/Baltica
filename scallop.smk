# -*- coding: utf-8
"""
Created on 17:07 27/07/2018 
Snakemake file for scalopadvising.


"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2018, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

names, = glob_wildcards('mappings/{names}.bam')

configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()
FILES = json.load(open(config['SAMPLES_JSON']))


rule all:
    input: 
        expand('scallop/{names}_denovo.gtf', names=names),
        'salmon/union.fa', 
        'salmon/salmon.index',   
        'salmon/quant.tsv.gz'


# rule optmise_scallop_parameters:
#     input: 'mappings/{names}.bam'
#     output: 'scallopadvising/{names}.gtf'
#     log:
#         'logs/{names}_scallopadvising.log'
#     params:
#         bin = '/home/tbrittoborges/bin/scallopadvising/ScallopAdvising.pl',
#         tmp = '/scratch/tbrittoborges/tmp/',
#         ref = '/biodb/genomes/homo_sapiens/GRCh38_90/GRCh38.90.gtf',
#         library_type = 'first'
#     shell:
#         """
#         module load scallop;
#         perl {params.bin} --log_file {log} --working_dir {params.tmp} --input_bam {input} --output_gtf {output} --reference {params.ref} --gtfcuff_path /home/tbrittoborges/bin/rnaseqtools/bin/gtfcuff --gffcompare_path /home/tbrittoborges/bin/gffcompare/gffcompare --library_type {params.library_type} /home/tbrittoborges/bin/scallopadvising/scallop_configs/*
#         """


rule scallop:
    input: 'mappings/{names}.bam'
    output: 'scallop/{names}.gtf'
    params:
        library_type = 'first'
    shell:
        """
        module load scallop;
        scallop -i {input} -o {output}  --library_type {params.library_type} --min_flank_length 6
        """


rule compare_to_reference:
    input: 'scallop/{names}.gtf'
    output: 'scallop/{names}_ref.gtf.tmap'
    params: 
        ref = '/biodb/genomes/homo_sapiens/GRCh38_90/GRCh38.90.gtf'
    shell: 'gffcompare -o {output} -r {params.ref} {input}'



rule extract_tx:
    params:
        ref_fa = '/prj/Niels_Gehring/newData_March_2018/tbb_analysis/Jenny_casc3/Homo_sapiens.GRCh38.dna_rm.toplevel.fa',
        ref = '/biodb/genomes/homo_sapiens/GRCh38_90/GRCh38.90.gtf'
    output: 'salmon/GRCh38_tx.fa'
    shell:
        """
        module load gffread
        gffread {params.ref} -g {params.ref_fa} -w {output}
        """


rule extract_denovo:
    input: 
        tmap = 'scallop/{names}_ref.gtf.tmap', 
        gtf = 'scallop/{names}.gtf'
    output: 'scallop/{names}_denovo.gtf'
    params: 
        ref = '/biodb/genomes/homo_sapiens/GRCh38_90/GRCh38.90.gtf'
    shell: """/home/tbrittoborges/bin/rnaseqtools/bin/gtfcuff puniq \
     {input.tmap} {input.gtf} {params.ref} {output}"""


# this step also concatenates all denovo gtf to one so we can a unique index for salmon
rule denovo_sequences:
    input:  expand('scallop/{names}_denovo.gtf', names=names)
    output: 'salmon/unique.fa'
    params:
        fasta = '/prj/Niels_Gehring/newData_March_2018/tbb_analysis/Jenny_casc3/Homo_sapiens.GRCh38.dna_rm.toplevel.fa'
    shell:  
        """
        module load gffread
        cat {input} > salmon/all.gff
        gffread salmon/all.gff -g {params.fasta} -w {output}
        """

rule merge_tx:
    input: 
        denovo_fa = 'salmon/unique.fa',
        ref_fa = 'salmon/GRCh38_tx.fa'
    output: 'salmon/union.fa'
    shell: "cat {input.denovo_fa} {input.ref_fa} > {output}"


rule salmon_index:
    input: 'salmon/union.fa'
    output: 'salmon/salmon.index'
    threads: 10
    log:
        'logs/salmon_index.log'
    shell : """
        module load salmon;
        salmon index -t {input} -i {output} -p {threads} &> {log}
        """


rule salmon_quant:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2'],
        index = 'salmon/salmon.index'
    output:
        'salmon/{sample}/quant.sf',
        'salmon/{sample}/lib_format_counts.json'
    log:
        'logs/{sample}_salmons_quant.log'
    shell: 
        """
        module load salmon;
        salmon quant -p 10 -i {input.index} \
         --libType A \
         -1 <(gunzip -c {input.r1}) \
         -2 <(gunzip -c {input.r2}) \
         -o salmon/{wildcards.sample}/ &> {log}
        """


rule collate_salmon:
    input:
        expand('salmon/{sample}/quant.sf', sample=names)
    output:
        'salmon/quant.tsv.gz'
    run:
        import gzip
        from os.path import basename, dirname

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = open(input[0]).readline()
            out.write(b('sample\t' + header))

            for i in input:
                sample = basename(dirname(i))
                lines = open(i)
                # Skip the header in each file.
                lines.readline()
                for line in lines:
                    out.write(b(sample + '\t' + line))