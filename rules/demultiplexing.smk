import os

rule demultiplex:
    output:
        temp(expand(os.path.join(config["general"]["output_dir"],"demultiplexed/{unit.sample}_{unit.unit}_R{read}.fastq.gz"), unit=units.reset_index().itertuples(), read=reads))
    params:
        filename = config["general"]["filename"].rstrip("/"),
        primertable = config["general"]["primertable"],
        demultiplexing = config["general"]["demultiplexing"],
        read_sorting = config['general']['read_sorting'],
        assembled = config['general']['already_assembled'],
        name_ext = config['merge']['name_ext'],
        output_dir = config['general']['output_dir']
    conda:
        "../envs/demultiplexing.yaml"
    script:
        "../scripts/demultiplexing.py"

rule unzip:
    input:
         os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_R{read}.fastq.gz")
    output:
        temp(os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{read}.tmp"))
    shell:
         "gunzip -c {input} > {output}"

rule check_format:
    input:
        os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{read}.tmp")
    output:
        temp(os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{read}.fastq"))
    shell:
        """
            if find {input} -not -type d -exec file '{{}}' ';' | grep CRLF
            then
                sed 's/\r$//' {input} > {output}
            else
                mv {input} {output}
            fi
        """
