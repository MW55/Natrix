rule fastqc:
    input:
        "demultiplexed/{sample}_{unit}_R{read}.fastq.gz"
    output:
        "results/qc/{sample}_{unit}_R{read}_fastqc.html",
        "results/qc/{sample}_{unit}_R{read}_fastqc.zip"
    threads: 15
    conda:
        "../envs/quality_control.yaml"
    shell:
        "fastqc -o results/qc/ {input} -t {threads}"

rule multiqc:
    input:
        expand("results/qc/{unit.sample}_{unit.unit}_R{read}_fastqc.zip",
        unit=units.reset_index().itertuples(), read=reads)
    output:
        "results/qc/multiqc_report.html"
    conda:
        "../envs/quality_control.yaml"
    shell:
        "multiqc results/qc/ -o results/qc/"