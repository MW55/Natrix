if config['general']['multiqc']:
    rule fastqc:
        input:
            "demultiplexed/{sample}_{unit}_R{read}.fastq.gz"
        output:
            "results/qc/{sample}_{unit}_R{read}_fastqc.html",
            "results/qc/{sample}_{unit}_R{read}_fastqc.zip"
        threads:
            config['general']['cores']
        conda:
            "../envs/quality_control.yaml"
        shell:
            "fastqc -o results/qc/ {input} -t {threads}"

    rule multiqc:
        input:
            expand("results/qc/{unit.sample}_{unit.unit}_R{read}_fastqc.zip",
            unit = units.reset_index().itertuples(), read=reads)
        output:
            "results/qc/multiqc_report.html"
        conda:
            "../envs/quality_control.yaml"
        shell:
            "multiqc results/qc/ -o results/qc/"

#I dont think this is really neccessary anymore, have to doublecheck
if config['general']['mid_check']:
    rule mid_check:
        input:
            data_folder = 'dereplicated',
            primer_table = config['general']['filename'] + '.csv'
        output:
            "logs/qc_done"
        log:
            "logs/quality_control.log"
        script:
            "../scripts/prim_mid_check.R"
