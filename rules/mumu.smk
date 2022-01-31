import os

rule generate_otu_fasta:
    input:
        expand(os.path.join(config["general"]["output_dir"], "finalData/{database}/swarm_mothur.csv"), database=config['classify']['database'])
    output:
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_mumu.fasta"), database=config['classify']['database'])
    script:
        "../scripts/generate_fasta.py"


rule vsearch_otu:
    input:
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_mumu.fasta"), database=config['classify']['database'])
    output:
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/match_scores.txt") , database=config['classify']['database'])
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
        "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"


rule run_mumu:
    input:
        os.path.join(config["general"]["output_dir"],"finalData/swarm_table.csv"),
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/match_scores.txt"), database=config['classify']['database'])
    output:
    	temp(expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.tmp"), database=config['classify']['database'])),
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.csv"), database=config['classify']['database'])
    log:
        os.path.join(config["general"]["output_dir"],"logs/finalData/otu_mumu.log")
    conda:
        "../envs/mumu.yaml"
    shell:
        """
        	cut -d "," -f 1,3- {input[0]} --output-delimiter="\t" > {output[0]};
        	mumu --otu_table {output[0]} --match_list {input[1]} --new_otu_table {output[1]} --log {log}
        """

rule merge_mumu_mothur_output:
    input:
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.csv"), database=config['classify']['database']),
        expand(os.path.join(config["general"]["output_dir"], "finalData/{database}/swarm_mothur.csv"), database=config['classify']['database'])
    output:
    	expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/FINAL_OUTPUT_OTU.csv"), database=config['classify']['database']),
    	expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/FINAL_OUTPUT_OTU_TABLE.csv"), database=config['classify']['database']),
    	expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/FINAL_OUTPUT_OTU_METADATA.csv"), database=config['classify']['database'])

    script:
            "../scripts/merge_mumu_output.py"
