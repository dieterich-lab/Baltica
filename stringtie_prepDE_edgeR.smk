# -*- coding: utf-8
"""
snakemake pipeline for stringtie to edgeR analysis

"""
configfile: "config.yml"
NAMES = config["samples"].keys()
SAMPLES = config["samples"].values()

rule all:
    input:
        'DEanalysis/gene_counts.csv'

rule create_sample_list:
    output:
        'DEanalysis/samples.txt'
    run:
        with open(str(output), 'w') as fou:
            for k, v in zip(NAMES, SAMPLES):
                fou.write('{} {}\n'.format(
                k, v.replace('/Aligned.noS.bam',
                '_StringTieBallgown/stringtieB_outfile.gtf')))

rule prepDE:
    input:
        'DEanalysis/samples.txt'
    output:
        g = 'DEanalysis/gene_counts.csv',
        t = 'DEanalysis/transcript_counts.csv'
    params:
        prepDE_path = '~/bin/prepDE.py',
        length = 75
    shell:
        '''
        module load python2
        python {params.prepDE_path} -i {input} --length={params.length} \
        -g {output.g} -t {output.t}
        '''

rule DE:
    input:
        rules.prepDE.output.g
    output:
        'DEanalysis/limmaReport/glMDSPlot/MDS-Plot.html',
        'DEanalysis/limmaReport/glMDPlot/MD-Plot.html'

    shell:
    '''
    module load R/3.4.1
    Rscript DE_analysis.R
    '''
