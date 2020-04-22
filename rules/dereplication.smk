rule cdhit:
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.fasta"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta",
        temp("results/assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta.clstr")
    conda:
        "../envs/dereplication.yaml"
    threads: config["general"]["cores"]
    params:
        id_percent=config["derep"]["clustering"],
        length_cutoff=config["derep"]["length_overlap"]
    shell:
        "cd-hit-est -i {input} -o {output[0]} -c {params.id_percent} -T"
        " {threads} -s {params.length_cutoff} -M 2000 -sc 1 -d 0"

rule cluster_sorting:
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta",
        "results/assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta.clstr",
        "results/assembly/{sample}_{unit}/{sample}_{unit}.fasta"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta"
    conda:
        "../envs/dereplication.yaml"
    params:
        repr=config["derep"]["representative"],
        length_cutoff=config["derep"]["length_overlap"]
    script:
        "../scripts/dereplication.py"
