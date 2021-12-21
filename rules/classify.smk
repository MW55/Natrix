import os

if config["classify"]["database"] == "pr2":
        rule mothur_classify:
            output:
                    expand(os.path.join(config["general"]["output_dir"], "results/finalData/pr2/mothur_out.summary")),
                    expand(os.path.join(config["general"]["output_dir"],"results/finalData/pr2/mothur_out.taxonomy"))
            input:
                    expand(os.path.join(config["general"]["output_dir"],"results/finalData/representatives.fasta")),
                    "database/pr2db.4.14.0.fasta"
            params:
                    template=config['database_path']['pr2_ref'],
                    taxonomy=config['database_path']['pr2_tax'],
                    search=config['classify']['search'],
                    method=config['classify']["method"],
                    threads=config['general']['cores'],
                    output=config['general']['output_dir']
            conda:
                    "../envs/mothur.yaml"
            log:
                    "results/logs/finalData/mothur_classify.log"
            shell:
                    """
                    mothur "#classify.seqs(fasta={input[0]}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})"; 
                    mv {params.output}/results/finalData/*.taxonomy {params.output}/results/finalData/pr2/mothur_out.taxonomy;
                    mv {params.output}/results/finalData/*.summary {params.output}/results/finalData/pr2/mothur_out.summary;
                    """

        rule mothur_swarm_merge:
            input:
                    swarm_abundance=os.path.join(config["general"]["output_dir"],"results/finalData/swarm_table.csv"),
                    mothur_tax=os.path.join(config["general"]["output_dir"],"results/finalData/pr2/mothur_out.taxonomy")
            output:
                    os.path.join(config["general"]["output_dir"],"results/finalData/pr2/swarm_mothur.csv")
            params:
                    scripts=config['scripts']['mothur_merge']
            shell:
                    "python {params.scripts} {input[0]} {input[1]} {output[0]}"

elif config["classify"]["database"] == "unite":
    rule mothur_classify:
        output:
            expand(os.path.join(config["general"]["output_dir"],"results/finalData/unite/mothur_out.summary")),
            expand(os.path.join(config["general"]["output_dir"],"results/finalData/unite/mothur_out.taxonomy"))
        input:
            expand(os.path.join(config["general"]["output_dir"],"results/finalData/representatives.fasta")),
            "database/unite_v8.3.fasta"
        params:
            template=config['database_path']['unite_ref'],
            taxonomy=config['database_path']['unite_tax'],
            search=config['classify']['search'],
            method=config['classify']["method"],
            output=config['general']['output_dir'],
            threads=config['general']['cores']
        conda:
            "../envs/mothur.yaml"
        log:
            "results/logs/finalData/mothur_classify.log"
        shell:
            """
            mothur "#classify.seqs(fasta={input[0]}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})"; 
            mv {params.output}/results/finalData/*.taxonomy {params.output}/results/finalData/unite/mothur_out.taxonomy;
            mv {params.output}/results/finalData/*.summary {params.output}/results/finalData/unite/mothur_out.summary;
            """

    rule mothur_swarm_merge:
        input:
            swarm_abundance=os.path.join(config["general"]["output_dir"],"results/finalData/swarm_table.csv"),
            mothur_tax=os.path.join(config["general"]["output_dir"],"results/finalData/unite/mothur_out.taxonomy")
        output:
            os.path.join(config["general"]["output_dir"],"results/finalData/unite/swarm_mothur.csv")
        params:
            scripts=config['scripts']['mothur_merge']
        shell:
            "python {params.scripts} {input[0]} {input[1]} {output[0]}"

elif config["classify"]["database"] == "silva":
    rule mothur_classify:
        output:
            expand(os.path.join(config["general"]["output_dir"],"results/finalData/silva/mothur_out.summary")),
            expand(os.path.join(config["general"]["output_dir"],"results/finalData/silva/mothur_out.taxonomy"))
        input:
            expand(os.path.join(config["general"]["output_dir"],"results/finalData/representatives.fasta")),
	    "database/silva_db.138.1.fasta", "database/silva_db.138.1.tax"
        params:
            template=config['database_path']['silva_ref'],
            taxonomy=config['database_path']['silva_tax'],
            search=config['classify']['search'],
            method=config['classify']["method"],
            output=config['general']['output_dir'],
            threads=config['general']['cores']
        conda:
            "../envs/mothur.yaml"
        log:
            os.path.join(config["general"]["output_dir"],"results/logs/finalData/mothur_classify.log")
        shell:
            """
            mothur "#classify.seqs(fasta={input[0]}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})"; 
            mv {params.output}/results/finalData/*.taxonomy {params.output}/results/finalData/silva/mothur_out.taxonomy;
            mv {params.output}/results/finalData/*.summary {params.output}/results/finalData/silva/mothur_out.summary;
            """

    rule mothur_swarm_merge:
        input:
            swarm_abundance=os.path.join(config["general"]["output_dir"],"results/finalData/swarm_table.csv"),
            mothur_tax=os.path.join(config["general"]["output_dir"],"results/finalData/silva/mothur_out.taxonomy")
        output:
            os.path.join(config["general"]["output_dir"],"results/finalData/silva/swarm_mothur.csv")
        params:
            scripts=config['scripts']['mothur_merge']
        shell:
            "python {params.scripts} {input[0]} {input[1]} {output[0]}"

