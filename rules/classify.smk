if config["classify"]["database"] == "PR2":
        rule mothur_classify:
            output:
                    expand("results/finalData/representatives.0.wang.tax.summary"), expand("results/finalData/representatives.1.wang.taxonomy")
            input:
                    expand("results/finalData/representatives.fasta")
            params:
                    template=config['database_path']['pr2_ref'],
                    taxonomy=config['database_path']['pr2_tax'],
                    search=config['classify']['search'],
                    method=config['classify']["method"],
                    threads=config['general']['cores']
            conda:
                    "../envs/mothur.yaml"
            log:
                    "results/logs/finalData/mothur_classify.log"
            shell:
                    """
                    mothur "#classify.seqs(fasta={input[0]}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})"; 
                    """

        rule mothur_swarm_merge:
            input:
                    swarm_abundance="results/finalData/swarm_table.csv",
                    mothur_tax="results/finalData/representatives.1.wang.taxonomy"
            output:
                    "results/finalData/swarm_mothur.csv"
            params:
                    scripts=config['scripts']['mothur_merge']
            shell:
                    "python {params.scripts} {input[0]} {input[1]} {output[0]}"

elif config["classify"]["database"] == "UNITE":
    rule mothur_classify:
        output:
            expand("results/finalData/representatives.0.wang.tax.summary"),expand("results/finalData/representatives.1.wang.taxonomy")
        input:
            expand("results/finalData/representatives.fasta")
        params:
            template=config['database_path']['unite_ref'],
            taxonomy=config['database_path']['unite_tax'],
            search=config['classify']['search'],
            method=config['classify']["method"],
            threads=config['general']['cores']
        conda:
            "../envs/mothur.yaml"
        log:
            "results/logs/finalData/mothur_classify.log"
        shell:
            """
            mothur "#classify.seqs(fasta={input[0]}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})"; 
            """

    rule mothur_swarm_merge:
        input:
            swarm_abundance="results/finalData/swarm_table.csv",
            mothur_tax="results/finalData/representatives.1.wang.taxonomy"
        output:
            "results/finalData/swarm_mothur.csv"
        params:
            scripts=config['scripts']['mothur_merge']
        shell:
            "python {params.scripts} {input[0]} {input[1]} {output[0]}"

elif config["classify"]["database"] == "SILVA":
    rule mothur_classify:
        output:
            expand("results/finalData/representatives.0.wang.tax.summary"),expand("results/finalData/representatives.1.wang.taxonomy")
        input:
            expand("results/finalData/representatives.fasta")
        params:
            template=config['database_path']['silva_ref'],
            taxonomy=config['database_path']['silva_tax'],
            search=config['classify']['search'],
            method=config['classify']["method"],
            threads=config['general']['cores']
        conda:
            "../envs/mothur.yaml"
        log:
            "results/logs/finalData/mothur_classify.log"
        shell:
            """
            mothur "#classify.seqs(fasta={input[0]}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})"; 
            """

    rule mothur_swarm_merge:
        input:
            swarm_abundance="results/finalData/swarm_table.csv",
            mothur_tax="results/finalData/representatives.1.wang.taxonomy"
        output:
            "results/finalData/swarm_mothur.csv"
        params:
            scripts=config['scripts']['mothur_merge']
        shell:
            "python {params.scripts} {input[0]} {input[1]} {output[0]}"

