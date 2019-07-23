rule all:
    input:
        "../database/ncbi/tax_lineage.h5"

rule create_blast_taxonomy:
    input:
        "../database/ncbi/fullnamelineage.dmp"
    output:
        "../database/ncbi/tax_lineage.h5"
        conda:
            "../envs/blast.yaml"
        log:
            "results/logs/finalData/BLAST.log"
        script:
            "../scripts/create_blast_taxonomy.py"
