rule demultiplex:
    input:
        samples=expand("input/{sample}_{unit}_R{read}.fastq.gz", sample=SAMPLES,unit=UNITS,read=READS),
        primertable="primertable.csv.tmp"
    output:
        temp(expand("demultiplexed/{sample}_{unit}_R{read}.fastq.gz",sample=SAMPLES,unit=UNITS,read=READS))
    params:
        demultiplexing = config["general"]["demultiplexing"],
        read_sorting = config['general']['read_sorting'],
        assembled = config['general']['already_assembled'],
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
