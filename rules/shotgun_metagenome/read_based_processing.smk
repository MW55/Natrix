# Taxonomic classification of shotgun metagenomic reads using Kraken2
# Depending on the experimental setting/question/resources, kraken supports different databases
# And overview and download links can be found at https://benlangmead.github.io/aws-indexes/k2

rule download_kraken2_bracken_db:
    output:
        directory("database/kraken2/"),
        "database/kraken2/hash.k2d",
        *["database/kraken2/database{}mers.kmer_distrib".format(l) for l in [50, 75, 100, 150, 200, 250, 300]]
    params:
        db_url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20241228.tar.gz"
    shell:
        """
        mkdir -p database/kraken2
        if [ ! -f database/kraken2/hash.k2d ]; then
            wget -O database/kraken2/kraken2_db.tar.gz {params.db_url}
            tar -xzf database/kraken2/kraken2_db.tar.gz -C database/kraken2
        fi
        """

rule kraken2_classification:
    input:
        fq1="results/filtered/{sample}_{unit}_clean.1",
        fq2="results/filtered/{sample}_{unit}_clean.2",
        db="database/kraken2/hash.k2d"
    output:
        report="results/taxonomy/{sample}_{unit}_kraken2_report.txt",
        classified="results/taxonomy/{sample}_{unit}_classified_reads.txt",
        unclassified_1="results/taxonomy/{sample}_{unit}_unclassified_reads_1.txt",
        unclassified_2="results/taxonomy/{sample}_{unit}_unclassified_reads_2.txt"
    params:
        threads=config["shotgun"]["kraken2"]["threads"],
        options=config["shotgun"]["kraken2"].get("options", "--confidence 0.1")
    conda:
        "../../envs/kraken2.yaml"
    log:
        "results/logs/{sample}_{unit}/{sample}_{unit}_kraken2.log"
    shell:
        """
        kraken2 --db database/kraken2 \
                --threads {params.threads} \
                {params.options} \
                --memory-mapping \
                --use-names \
                --paired {input.fq1} {input.fq2} \
                --report {output.report} \
                --report-zero-counts \
                --output {output.classified} \
                --unclassified-out results/taxonomy/{wildcards.sample}_{wildcards.unit}_unclassified_reads#.txt \
                2> {log}
        """

# REMOVE memory-mapping if not running locally! It reads from disk, not memory

rule bracken_abundance:
    input:
        kraken_report="results/taxonomy/{sample}_{unit}_kraken2_report.txt",
        kmer_distrib=lambda wildcards: f"database/kraken2/database{config['shotgun']['bracken']['read_length']}mers.kmer_distrib"  # Use correct read length
    output:
        bracken_report="results/taxonomy/{sample}_{unit}_bracken_report.txt"
    params:
        db="database/kraken2",
        threshold=config["shotgun"]["bracken"]["threshold"],
        threads=config["shotgun"]["bracken"]["threads"],
        read_length=config["shotgun"]["bracken"]["read_length"],
        level=config["shotgun"]["bracken"]["level"]
    conda:
        "../../envs/kraken2.yaml"
    log:
        "results/logs/{sample}_{unit}/{sample}_{unit}_bracken.log"
    shell:
        """
        bracken -d {params.db} \
                -i {input.kraken_report} \
                -o {output.bracken_report} \
                -r {params.read_length} \
                -l {params.level} \
                -t {params.threshold} \
                -w {params.threads} 2> {log}
        """
        #-l is the level, C is class, i.e. very high level, S is species. I used C for the example_data, as the classification did not result in much