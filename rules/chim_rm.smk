import os

rule vsearch:
    input:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta") if config["general"]["seq_rep"] == "OTU" else os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_dada.fasta")
    output:
        uchime_out=temp(os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.vsearch.txt")),
        chim=os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.chimera.fasta"),
        nonchim=os.path.join(config["general"]["output_dir"],"finalData/{sample}_{unit}.nonchimera.fasta")
    params:
        beta=config["chim"]["beta"],
        pseudo_c=config["chim"]["pseudo_count"],
        abskew=config["chim"]["abskew"]
    threads: config["general"]["cores"]
    conda:
         "../envs/vsearch.yaml"
    log:
        os.path.join(config["general"]["output_dir"],"logs/{sample}_{unit}/vsearch.log")
    shell:
        "vsearch --uchime3_denovo {input} -uchimeout {output.uchime_out}"
        " -chimeras {output.chim} -nonchimeras {output.nonchim} -xn {params.beta}"
        "  -dn {params.pseudo_c} -abskew {params.abskew} --log {log} 2>&1"
