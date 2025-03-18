rule download_card:
    output:
        card_archive="database/CARD/data"
    shell:
        """
        mkdir -p database/CARD
        if [ ! -f {output.card_archive} ]; then
            wget -O {output.card_archive} https://card.mcmaster.ca/latest/data
        fi
        """

rule extract_card:
    input:
        card_archive="database/CARD/data"
    output:
        card_json="database/CARD/card.json"
    shell:
        """
        tar -xvf {input.card_archive} -C database/CARD ./card.json
        """

rule create_annotation:
    input:
        card_json="database/CARD/card.json",
    output:
        card_db_fasta="database/CARD/card_database_v4.0.0.fasta",
    conda:
        "../../envs/rgi.yaml"
    shell:
        """
        rgi card_annotation -i {input.card_json}
        mv card_database_v4.0.0.fasta {output.card_db_fasta}
        """

rule load_card_db:
    input:
        card_json = "database/CARD/card.json",
        card_db_fasta="database/CARD/card_database_v4.0.0.fasta",
    output:
        touch("database/CARD/db_loaded.done")
    conda:
        "../../envs/rgi.yaml"
    shell:
        """
        rgi load --card_json {input.card_json} --card_annotation {input.card_db_fasta} --local
        touch {output}
        """

rule rgi:
    input:
        r1="results/filtered/{sample}_{unit}_clean.1",
        r2=lambda wildcards: "results/filtered/{sample}_{unit}_clean.2" if config["merge"]["paired_End"] else [],
        db_done="database/CARD/db_loaded.done"
    output:
        "results/rgi/{sample}_{unit}.overall_mapping_stats.txt"
    params:
        threads=2
    conda:
        "../../envs/rgi.yaml"
    log:
        "results/logs/{sample}_{unit}/rgi.log"
    script:
        "../../scripts/resistance_prediction/rgi.py"

rule merge_rgi_reports:
    input:
        reports = expand("results/rgi/{unit.sample}_{unit.unit}.overall_mapping_stats.txt", unit=units.reset_index().itertuples())
    output:
        "results/rgi/merged_rgi.txt"
    log:
        "results/logs/merge_rgi.log"
    shell:
        "cat {input.reports} > {output}"