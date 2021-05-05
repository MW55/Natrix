rule fastqc:
    input:
        "demultiplexed/{sample}_{unit}_R{read}.fastq.gz"
    output:
        "results/qc/{sample}_{unit}_R{read}_fastqc.html",
        "results/qc/{sample}_{unit}_R{read}_fastqc.zip"
    threads: 20
    conda:
        "../envs/quality_control.yaml"
    shell:
        "fastqc -o results/qc/ {input} -t {threads}"

rule multiqc:
    input:
        expand("results/qc/{{sample}}_{{unit}}_R{read}_fastqc.zip",read=READS)
    output:
        "results/qc/multiqc_report.html"
    conda:
        "../envs/quality_control.yaml"
    shell:
        "multiqc results/qc/ -o results/qc/"