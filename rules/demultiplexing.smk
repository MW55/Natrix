rule demultiplex:
    output:
        (expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_assembled.fastq", unit=units.reset_index().itertuples()) 
        if config["general"]["already_assembled"] 
        else temp(expand("demultiplexed/{unit.sample}_{unit.unit}_R{read}.fastq.gz", unit=units.reset_index().itertuples(), read=reads)) )
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
