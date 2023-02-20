import os

if config['classify']['mothur']:
    rule generate_otu_fasta:
        input:
            os.path.join(config["general"]["output_dir"], "finalData/{database}/full_table.csv")
        output:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_mumu.fasta")
        script:
            "../scripts/generate_fasta.py"

    rule vsearch_otu:
        input:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_mumu.fasta")
        output:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/match_scores.txt")
        conda:
            "../envs/vsearch.yaml"
        shell:
            "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
            "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"

    rule run_mumu:
        input:
            os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" and config["dataset"]["nanopore"]== "FALSE" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
            expand(os.path.join(config["general"]["output_dir"],"mothur/{database}/match_scores.txt"), database=config['classify']['database']),
        output:
            temp(expand(os.path.join(config["general"]["output_dir"], "mothur/{database}/OTU_table_mumu.tmp"), database=config['classify']['database'])),
            expand(os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_table_mumu.csv"), database=config['classify']['database'])
        log:
            os.path.join(config["general"]["output_dir"],"logs/otu_mumu.log")
        conda:
            "../envs/mumu.yaml"
        shell:
            """
            	cut -d "," -f 1,3- {input[0]} --output-delimiter="\t" > {output[0]};
            	mumu --otu_table {output[0]} --match_list {input[1]} --new_otu_table {output[1]} --log {log}
            """

    rule merge_mumu_mothur_output:
        input:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_table_mumu.csv"),
        	os.path.join(config["general"]["output_dir"],"finalData/{database}/full_table.csv"),
        output:
            os.path.join(config["general"]["output_dir"],"finalData/{database}/full_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/metadata_table_mumu.csv")
        script:
            "../scripts/merge_mumu_output.py"

else:
    rule generate_otu_fasta:
        input:
            expand(os.path.join(config["general"]["output_dir"], "finalData/blast_{database}/full_table.csv"), database=config['blast']['database'].lower())
        output:
            os.path.join(config["general"]["output_dir"],"blast/OTU_mumu.fasta")
        script:
            "../scripts/generate_fasta.py"

    rule vsearch_otu:
        input:
            os.path.join(config["general"]["output_dir"],"blast/OTU_mumu.fasta")
        output:
            os.path.join(config["general"]["output_dir"],"blast/match_scores.txt")
        conda:
            "../envs/vsearch.yaml"
        shell:
            "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
            "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"

    rule run_mumu:
        input:
            os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" and config["dataset"]["nanopore"]== "FALSE" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
            expand(os.path.join(config["general"]["output_dir"],"blast/match_scores.txt"), database=config['classify']['database']),
        output:
            temp(expand(os.path.join(config["general"]["output_dir"], "blast/OTU_table_mumu.tmp"), database=config['classify']['database'])),
            expand(os.path.join(config["general"]["output_dir"],"blast/OTU_table_mumu.csv"), database=config['classify']['database'])
        log:
            os.path.join(config["general"]["output_dir"],"logs/otu_mumu.log")
        conda:
            "../envs/mumu.yaml"
        shell:
            """
            	cut -d "," -f 1,3- {input[0]} --output-delimiter="\t" > {output[0]};
            	mumu --otu_table {output[0]} --match_list {input[1]} --new_otu_table {output[1]} --log {log}
            """

    rule merge_mumu_blast_output:
        input:
            os.path.join(config["general"]["output_dir"],"blast/OTU_table_mumu.csv"),
        	os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/full_table.csv"),
        output:
            os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/full_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/OTU_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/metadata_table_mumu.csv")
        script:
            "../scripts/merge_mumu_output2.py"
