import os

if config["classify"]["database"] == "pr2":
        rule mothur_classify:
            output:
                    expand(os.path.join(config["general"]["output_dir"], "finalData/pr2/mothur_out.summary")),
                    expand(os.path.join(config["general"]["output_dir"],"finalData/pr2/mothur_out.taxonomy"))
            input:
                    expand(os.path.join(config["general"]["output_dir"],"finalData/representatives.fasta")),
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
                    os.path.join(config["general"]["output_dir"], "logs/finalData/mothur_classify.log")
            shell:
                    """
                        mothur "#set.logfile(name={log}); classify.seqs(fasta={input[0]}, cutoff=0, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                        mv {params.output}/finalData/*.taxonomy {params.output}/finalData/pr2/mothur_out.taxonomy;
                        mv {params.output}/finalData/*.summary {params.output}/finalData/pr2/mothur_out.summary;
                    """

        rule mothur_swarm_merge:
            input:
                    swarm_abundance=os.path.join(config["general"]["output_dir"],"finalData/swarm_table.csv"),
                    mothur_tax=os.path.join(config["general"]["output_dir"],"finalData/pr2/mothur_out.taxonomy")
            output:
                    os.path.join(config["general"]["output_dir"],"finalData/pr2/swarm_mothur.csv")
            script:
                   "../scripts/merge_mothur_with_swarm.py"


elif config["classify"]["database"] == "unite":
    rule mothur_classify:
        output:
            expand(os.path.join(config["general"]["output_dir"],"finalData/unite/mothur_out.summary")),
            expand(os.path.join(config["general"]["output_dir"],"finalData/unite/mothur_out.taxonomy"))
        input:
            expand(os.path.join(config["general"]["output_dir"],"finalData/representatives.fasta")),
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
            "logs/finalData/mothur_classify.log"
        shell:
            """
                mothur "#classify.seqs(fasta={input[0]}, cutoff=0, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                mv {params.output}/finalData/*.taxonomy {params.output}/finalData/unite/mothur_out.taxonomy;
                mv {params.output}/finalData/*.summary {params.output}/finalData/unite/mothur_out.summary;
            """

    rule mothur_swarm_merge:
        input:
            swarm_abundance=os.path.join(config["general"]["output_dir"],"finalData/swarm_table.csv"),
            mothur_tax=os.path.join(config["general"]["output_dir"],"finalData/unite/mothur_out.taxonomy")
        output:
            os.path.join(config["general"]["output_dir"],"finalData/unite/swarm_mothur.csv")
        script:
            "../scripts/merge_mothur_with_swarm.py"

elif config["classify"]["database"] == "silva":
    rule mothur_classify:
        output:
            expand(os.path.join(config["general"]["output_dir"],"finalData/silva/mothur_out.summary")),
            expand(os.path.join(config["general"]["output_dir"],"finalData/silva/mothur_out.taxonomy"))
        input:
            expand(os.path.join(config["general"]["output_dir"],"finalData/representatives.fasta")),
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
            os.path.join(config["general"]["output_dir"],"logs/finalData/mothur_classify.log")
        shell:
            """
                mothur "#classify.seqs(fasta={input[0]}, cutoff=0, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                mv {params.output}/finalData/*.taxonomy {params.output}/finalData/silva/mothur_out.taxonomy;
                mv {params.output}/finalData/*.summary {params.output}/finalData/silva/mothur_out.summary;
            """

    rule mothur_swarm_merge:
        input:
            swarm_abundance=os.path.join(config["general"]["output_dir"],"finalData/swarm_table.csv"),
            mothur_tax=os.path.join(config["general"]["output_dir"],"finalData/silva/mothur_out.taxonomy")
        output:
            os.path.join(config["general"]["output_dir"],"finalData/silva/swarm_mothur.csv")
        script:
            "../scripts/merge_mothur_with_swarm.py"
