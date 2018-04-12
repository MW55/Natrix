if config['general']['multiqc']:
    rule fastqc:
        input:
            data_folder + "/{sample}_{unit}_R{group}.fastq.gz"
        output:
            "results/qc/{sample}_{unit}_R{group}_fastqc.html",
            "results/qc/{sample}_{unit}_R{group}_fastqc.zip",
        shell:
            "fastqc -o results/qc/ {input}"

    rule multiqc:
        input:
            expand("results/qc/{unit.sample}_{unit.unit}_R{group}_fastqc.zip",
            unit = units.reset_index().itertuples(), group=groups)
        output:
            "results/qc/multiqc_report.html"
        shell:
            "multiqc results/qc/ -o results/qc/" 

if config['general']['mid_check']:
    rule mid_check:
        input:
            data_folder,
            primer_table = data_folder + '.csv'
        output:
            "logs/qc_done"
        log:
            "logs/quality_control.log"
        script:
            "../scripts/prim_mid_check.R"
