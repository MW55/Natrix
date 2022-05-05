import os

rule cdhit:
    input:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.fasta")
    output:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta"),
        temp(os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta.clstr"))
    conda:
        "../envs/dereplication.yaml"
    threads: config["general"]["cores"]
    params:
        id_percent=config["derep"]["clustering"],
        length_cutoff=config["derep"]["length_overlap"],
        memory=config["general"]["memory"]
    shell:
        "cd-hit-est -i {input} -o {output[0]} -c {params.id_percent} -T"
        " {threads} -s {params.length_cutoff} -M {params.memory} -sc 1 -d 0"

rule cluster_sorting:
    input:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta"),
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta.clstr"),
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.fasta")
    output:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta")
    conda:
        "../envs/dereplication.yaml"
    params:
        repr=config["derep"]["representative"],
        length_cutoff=config["derep"]["length_overlap"]
    script:
        "../scripts/dereplication.py"
