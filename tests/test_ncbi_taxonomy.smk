rule all:
    input:
        "files/blast_taxonomic_lineage.tsv"

rule ncbi_taxonomy:
    input:
        blast_result = "files/blast_taxonomy.tsv",
        lineage = "../database/ncbi/tax_lineage.h5"
    output:
        tax_lineage = "files/blast_taxonomic_lineage.tsv",
        all_tax = "files/blast_taxonomy_all.tsv"
    params:
        max_target_seqs=10
    conda:
        "../envs/blast.yaml"
    log:
        "results/logs/finalData/BLAST.log"
    script:
        "../scripts/ncbi_taxonomy.py"
