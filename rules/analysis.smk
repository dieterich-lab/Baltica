
workdir: config.get("path", ".")

rule parse_majiq:
    input:
    output:
    shell:
        "Rscript --vanilla analysis/parse_majiq_output.R"


rule parse_leafcutter:
    input:
    output:
    shell:
        "Rscript --vanilla analysis/parse_leafcutter_results.R leafcutter"

rule parse_majiq:
    input:
    output:
    shell:
        "Rscript analysis/parse_junctionseq_output.R"


rule annotate:
    input:
        table = "{method}/{method}_junctions.csv",
        annotatio = "merged/merged.combined.gtf"
    output:
        "{method}/{method}_junctions_annotated.csv"
    shell:
        "Rscript analysis/annotate_SJ.R -i {input.table} -a {input.annotation} -o {output}"Rscript analysis/annotate_SJ.R -i majiq/majiq_junctions.csv -a denovo_tx/merged/merged.combined.gtf -o majiq/majiq_junctions_annotated.csv

rule write_simplified_result:

    shell: Rscript Sch4_7_simplify_DJU.R
