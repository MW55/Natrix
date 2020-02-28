rule vsearch:
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta"
    output:
        uchime_out=temp("results/assembly/{sample}_{unit}/{sample}_{unit}.vsearch.txt"),
        chim="results/assembly/{sample}_{unit}/{sample}_{unit}.chimera.fasta",
        nonchim="results/finalData/{sample}_{unit}.nonchimera.fasta"
    params:
        beta=config["chim"]["beta"],
        pseudo_c=config["chim"]["pseudo_count"],
        abskew=config["chim"]["abskew"]
    threads: 20
    conda:
         "../envs/vsearch.yaml"                                                                  
    shell:
        "vsearch --uchime3_denovo {input} -uchimeout {output.uchime_out}"
        " -chimeras {output.chim} -nonchimeras {output.nonchim} -xn {params.beta}"
        "  -dn {params.pseudo_c} -abskew {params.abskew} 2>&1"
