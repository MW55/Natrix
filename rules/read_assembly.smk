def get_fastq(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        return expand("demultiplexed/{sample}_{unit}_{group}.fastq",
                        group=[1,2], **wildcards)
    return "demultiplexed/{sample}_{unit}_1.fastq".format(**wildcards)

rule demultiplex:
    #input:
         #expand(config["general"]["filename"] + "/{sample}_{unit}_R{read}.fastq.gz")
    output:
        expand("demultiplexed/{unit.sample}_{unit.unit}_R{read}.fastq.gz", unit=units.reset_index().itertuples(), read=reads)
    params:
        filename = config["general"]["filename"],
        primertable = config["general"]["primertable"],
        demultiplexing = config["general"]["demultiplexing"],
        read_sorting = config['general']['read_sorting'],
        assembled = config['general']['already_assembled'],
        name_ext = config['merge']['name_ext']
    conda:
        "../envs/demultiplexing.yaml"
    script:
        "../scripts/demultiplexing.py"

rule unzip:
    input:
         "demultiplexed/{sample}_{unit}_R{read}.fastq.gz"
    output:
        temp("demultiplexed/{sample}_{unit}_{read}.fastq")
    shell: 
         "gunzip -c {input} > {output}"

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

rule DADA2:
    input:
        forward = expand(
            "results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_1_cut.fastq", unit=units.reset_index().itertuples()),
        reverse = expand(
            "results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_2_cut.fastq", unit=units.reset_index().itertuples()) if config["merge"]["paired_End"] == True else []
    output:
        expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_dada.fasta", unit=units.reset_index().itertuples())
    params:
        paired_end=config["merge"]["paired_End"]
    conda:
        "../envs/dada2.yaml"
    log:
        "results/logs/finalData/DADA2.log"
    script:
        "../scripts/dada2.R"


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