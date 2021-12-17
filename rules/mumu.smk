import os

rule generate_otu_fasta:
    input:
        expand(os.path.join(config["general"]["output_dir"], "results/finalData/{database}/swarm_mothur.csv"), database=config['classify']['database'])
    output:
        expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/OTU_fasta.txt"), database=config['classify']['database'])
    script:
        "../scripts/generate_fasta.py"


rule vsearch_otu:
    input:
        expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/OTU_fasta.txt"), database=config['classify']['database'])
    output:
        expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/match_scores.txt") , database=config['classify']['database'])
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
        "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"


rule run_mumu:
    input:
        expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/swarm_mothur.csv"), database=config['classify']['database']),
        expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/match_scores.txt"), database=config['classify']['database'])
    output:
        expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/OTU_table_mumu.csv"), database=config['classify']['database'])
    log:
        os.path.join(config["general"]["output_dir"],"results/logs/finalData/otu_mumu.log")
    conda:
        "../envs/mumu.yaml"
    shell:
        "mumu --otu_table {input[0]} --match_list {input[1]} --new_otu_table {output} --log {log}"


rule edit_output:
     input:
           expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/OTU_table_mumu.csv"), database=config['classify']['database'])
     output:
           expand(os.path.join(config["general"]["output_dir"],"results/finalData/{database}/FINAL_OUTPUT_OTU.txt"), database=config['classify']['database'])
     script:
            "../scripts/edit_mumu_output.py"

