# Functional annotation of shotgun metagenomic reads using HUMAnN3.

rule download_humann3_db:
    output:
        directory("database/humann3/"),
        "database/humann3/chocophlan.tar.gz",
        "database/humann3/uniref90.tar.gz"
    params:
        chocophlan_url="http://cmprod1.cibio.unitn.it/databases/HUMAnN/full_chocophlan.v296_201901.tar.gz", #"https://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.tar.gz",
        uniref_url="http://cmprod1.cibio.unitn.it/databases/HUMAnN/uniref90_annotated_v201901.tar.gz" #"https://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_ec_filtered/uniref90_ec_filtered.tar.gz"
    shell:
        """
        mkdir -p database/humann3
        if [ ! -d database/humann3/chocophlan ]; then
            wget -O database/humann3/chocophlan.tar.gz {params.chocophlan_url}
            tar -xzf database/humann3/chocophlan.tar.gz -C database/humann3/
        fi
        if [ ! -d database/humann3/uniref90 ]; then
            wget -O database/humann3/uniref90.tar.gz {params.uniref_url}
            tar -xzf database/humann3/uniref90.tar.gz -C database/humann3/
        fi
        """

rule concatenate_reads:
    input:
        fq1="results/filtered/{sample}_{unit}_clean.1",
        fq2="results/filtered/{sample}_{unit}_clean.2"
    output:
        merged_fq="results/filtered/{sample}_{unit}_merged.fastq"
    shell:
        """
        cat {input.fq1} {input.fq2} > {output.merged_fq}
        """

rule humann3_functional_annotation:
    input:
        merged_fq="results/filtered/{sample}_{unit}_merged.fastq",
        taxonomy_report="results/taxonomy/{sample}_{unit}_kraken2_report.txt",
        chocophlan_db="database/humann3/chocophlan.tar.gz",
        uniref_db="database/humann3/uniref90.tar.gz"
    output:
        gene_families="results/functional/{sample}_{unit}_genefamilies.tsv",
        pathways="results/functional/{sample}_{unit}_pathways.tsv",
        log="results/logs/{sample}_{unit}/{sample}_{unit}_humann3.log"
    params:
        threads=config["shotgun"]["humann3"]["threads"],
        out_dir="results/functional/{sample}_{unit}"
    conda:
        "../../envs/humann3.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        humann --input {input.merged_fq} \
               --taxonomic-profile {input.taxonomy_report} \
               --output {params.out_dir} \
               --output-basename {wildcards.sample}_{wildcards.unit} \
               --nucleotide-database {input.chocophlan_db} \
               --protein-database {input.uniref_db} \
               --threads {params.threads} &> {output.log}
        """

rule humann3_renorm:
    input:
        gene_families="results/functional/{sample}_{unit}_genefamilies.tsv",
        pathways="results/functional/{sample}_{unit}_pathways.tsv"
    output:
        gene_families_relab="results/functional/{sample}_{unit}_genefamilies_relab.tsv",
        pathways_relab="results/functional/{sample}_{unit}_pathways_relab.tsv"
    conda:
        "../../envs/humann3.yaml"
    shell:
        """
        humann_renorm_table --input {input.gene_families} --output {output.gene_families_relab} --units relab
        humann_renorm_table --input {input.pathways} --output {output.pathways_relab} --units relab
        """

rule humann3_metacyc_pathways:
    input:
        pathways_relab="results/functional/{sample}_{unit}_pathways_relab.tsv"
    output:
        metacyc_pathways="results/functional/{sample}_{unit}_metacyc_pathways.tsv"
    conda:
        "../../envs/humann3.yaml"
    shell:
        """
        humann_barplot --input {input.pathways_relab} --output {output.metacyc_pathways} --level metacyc
        """