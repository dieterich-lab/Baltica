comp = config["contrasts"].keys()
workdir: config['path']


# snakemake --dag -n | dot -Tsvg > dag.svg

subworkflow qc:
    workdir:
        config['path']
    snakefile:
        "qc.smk"
    configfile: 
        "/home/tbrittoborges/data/config.yml"

subworkflow junctionseq:
    workdir:
        config['path']
    snakefile:
        "junctionseq.smk"
    configfile: 
        "/home/tbrittoborges/Baltica/data/config.yml"

subworkflow leafcutter:
    workdir:
        config['path']
    snakefile:
        "leafcutter.smk"
    configfile: 
        "/home/tbrittoborges/Baltica/data/config.yml"

subworkflow majiq:
    workdir:
        config['path']
    snakefile:
        "majiq.smk"
    configfile: 
        "/home/tbrittoborges/Baltica/data/config.yml"

subworkflow stringtie:
    workdir:
        config['path']
    snakefile:
        "stringtie.smk"
    configfile: 
        "/home/tbrittoborges/Baltica/data/config.yml"

subworkflow rmats:
    workdir:
        config['path']
    snakefile:
        "rmats.smk"
    configfile: 
        "/home/tbrittoborges/Baltica/data/config.yml"

subworkflow analysis:
    workdir:
        config['path']
    snakefile:
        "analysis.smk"
    configfile: 
        "/home/tbrittoborges/Baltica/data/config.yml"

rule final:
    input:
        majiq(
            expand("majiq/voila/{comp}_voila.tsv",comp=comp)),
        rmats(
            expand("rmats/{comp}/",comp=comp)),
        leafcutter(
            expand("leafcutter/{comp}/{comp}_cluster_significance.txt",comp=comp)),
        junctionseq(
            expand("junctionseq/analysis/{comp}_sigGenes.results.txt.gz", comp=comp)),
        stringtie("stringtie/merged/merged.combined.gtf"),
        analysis(
            "results/SJ_annotated_assigned_simple.xlsx"),    
