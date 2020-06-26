rule all:
    input: 'results/SJ_annotated_assigned.csv'

subworkflow qc:
    workdir:
        config.get("path", ".")
    snakefile:
        "rules/qc.smk"

subworkflow majiq:
    workdir:
        config.get("path", ".")
    snakefile:
        "rules/qc.smk"

subworkflow leafcutter:
    workdir:
        config.get("path", ".")
    snakefile:
        "rules/qc.smk"

subworkflow junctionseq:
    workdir:
        config.get("path", ".")
    snakefile:
        "rules/qc.smk"

subworkflow analysis:
    workdir:
        config.get("path", ".")
    snakefile:
        "rules/analysis.smk"

