rule all:
    input:
        "files/filtered_blast_table.csv"

rule merge_results:
    input:
        merged_swarm="files/swarm_table.csv",
        blast_result="files/blast_taxonomic_lineage.tsv"
    output:
        "files/filtered_blast_table.csv"
    conda:
        "../envs/merge_results.yaml"
    log:
        "results/logs/finalData/BLAST.log"
    script:
        "../scripts/merge_results.py"
