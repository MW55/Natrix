import os

rule fastqc:
    input:
        os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_R{read}.fastq.gz")
    output:
        os.path.join(config["general"]["output_dir"],"qc/{sample}_{unit}_R{read}_fastqc.html"),
        os.path.join(config["general"]["output_dir"],"qc/{sample}_{unit}_R{read}_fastqc.zip")
    params:
        config["general"]["output_dir"]
    threads: 20
    conda:
        "../envs/quality_control.yaml"
    shell:
        "fastqc -o {params}/qc/ {input} -t {threads}"

rule multiqc:
    input:
        expand(os.path.join(config["general"]["output_dir"],"qc/{unit.sample}_{unit.unit}_R{read}_fastqc.zip"),
        unit=units.reset_index().itertuples(), read=reads)
    output:
        os.path.join(config["general"]["output_dir"],"qc/multiqc_report.html")
    params:
        config["general"]["output_dir"]
    conda:
        "../envs/quality_control.yaml"
    shell:
        "multiqc {params}/qc/ -o {params}/qc/"
