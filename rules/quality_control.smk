rule fastqc:
    input:
        "demultiplexed/{sample}_{unit}_{unit.group}.fastq.gz"
    output:
        "results/qc/{sample}_{unit}_{unit.group}_fastqc.html",
        "results/qc/{sample}_{unit}_{unit.group}_fastqc.zip"
    threads:
        config['general']['cores']
    conda:
        "../envs/quality_control.yaml"
    shell:
        "fastqc -o results/qc/ {input} -t {threads}"

rule multiqc:
    input:
        'results/qc/'
    output:
        "results/qc/multiqc_report.html"
    conda:
        "../envs/quality_control.yaml"
    shell:
        "multiqc {input} -o {input}"
