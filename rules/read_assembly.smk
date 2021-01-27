def get_fastq(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        return expand("demultiplexed/{sample}_{unit}_{group}.fastq",
                        group=[1,2], **wildcards)
    return "demultiplexed/{sample}_{unit}_1.fastq".format(**wildcards)

rule define_primer:
    input:
        primer_table=config["general"]["primertable"]
    output:
        "primer_table.csv"
    params:
       paired_end=config["merge"]["paired_End"],
       offset=config["qc"]["primer_offset"],
       bar_removed=config["qc"]["barcode_removed"],
       all_removed=config["qc"]["all_primer"]
    conda:
        "../envs/define_primer.yaml"
    script:
        "../scripts/define_primer.py"

rule prinseq:
    input:
        sample=get_fastq
    output:
        expand(
        "results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq",
        read=reads)
    params:
        config["qc"]["mq"]
    log:
        "results/logs/{sample}_{unit}/prinseq.log"
    conda:
        "../envs/prinseq.yaml"
    script:
        "../scripts/prinseq.py"

rule cutadapt:
    input:
        expand(
        "results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq",
        read=reads),
        primer_t="primer_table.csv"
    output:
        expand(
        "results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}_cut.fastq",
        read=reads)
    params:
        paired_end=config["merge"]["paired_End"],
        bar_removed=config["qc"]["barcode_removed"],
        prim_rm=config["qc"]["all_primer"],
        minlen=config["qc"]["minlen"],
        maxlen=config["qc"]["maxlen"]
    conda:
        "../envs/cutadapt.yaml"
    log:
        "results/logs/{sample}_{unit}/cutadapt.log"
    script:
        "../scripts/cutadapt.py"

rule assembly:
    input:
        expand(
        "results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq",
        read=reads),
        primer_t="primer_table.csv"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq"
    threads: 20
    params:
        paired_end=config["merge"]["paired_End"],
        threshold=config["qc"]["threshold"],
        minoverlap=config["qc"]["minoverlap"],
        minlen=config["qc"]["minlen"],
        maxlen=config["qc"]["maxlen"],
        minqual=config["qc"]["minqual"],
        prim_rm=config["qc"]["all_primer"]
    conda:
        "../envs/assembly.yaml"
    log:
        "results/logs/{sample}_{unit}/read_assembly.log"
    script:
        "../scripts/assembly.py"

rule copy_to_fasta:
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.fasta"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -a {input} > {output}"
