import os

if config["blast"]["database"] == "SILVA":

    rule make_silva_db:
        output:
            expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".nhr", ".nin", ".nog", ".nsd", ".nsi", ".nsq"]),
            config["blast"]["db_path"] + ".fasta"
        params:
            db_path=config["blast"]["db_path"]
        conda:
            "../envs/blast.yaml"
        shell:
            """
                dir_name=$(dirname {params[0]});
                wget -P $dir_name/ --progress=bar https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_132_SSURef_tax_silva.fasta.gz;
                gunzip -c $dir_name/SILVA_132_SSURef_tax_silva.fasta.gz > $dir_name/silva.db.fasta;
                makeblastdb -in $dir_name/silva.db.fasta -dbtype nucl -parse_seqids -out $dir_name/silva.db
            """

    rule create_silva_taxonomy:
        input: config["blast"]["db_path"] + ".fasta"
        output: os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        conda:
            "../envs/blast.yaml"
        script:
            "../scripts/create_silva_taxonomy.py"

elif config["blast"]["database"] == "NCBI":

    rule make_ncbi_db:
        output:
            expand(config["blast"]["db_path"] + ".00" + "{file_extension}", file_extension=[".nhr", ".nin", ".nog", ".nsd", ".nsi", ".nsq"]),
            listing = temp(".listing"),
            path = config["blast"]["db_path"]
        params:
            db_path=config["blast"]["db_path"]
        conda:
            "../envs/blast.yaml"
        shell:
            """
                echo 'what--------------';
                dir_name=$(dirname {params[0]});
                wget --spider --no-remove-listing ftp://ftp.ncbi.nlm.nih.gov/blast/db/
                number=$(awk '$9 ~ /^nt.[0-9]*.tar.gz[^.]/ {{print substr($9,4,2)}}' {output.listing} | tail -n 1)
                for i in `seq -w 00 $number`;
                    do wget -N -P $dir_name/ --progress=bar ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.${{i}}.tar.gz;
                tar xvzf $dir_name/nt.${{i}}.tar.gz -C $dir_name;
                done;
                touch {output.path}
            """

    rule download_taxonomy:
        output:
            os.path.join(os.path.dirname(config["blast"]["db_path"]), "fullnamelineage.dmp")
        params:
            db_path=config["blast"]["db_path"]
        conda:
            "../envs/blast.yaml"
        shell:
            """
                echo 'what--------------';
                dir_name=$(dirname {params[0]});
                wget -N -P $dir_name/ --progress=bar ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz;
                wget -N -P $dir_name/ --progress=bar ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz;
                tar xvzf $dir_name/new_taxdump.tar.gz -C $dir_name;
            """

    rule create_blast_taxonomy:
        input:
            os.path.join(os.path.dirname(config["blast"]["db_path"]), "fullnamelineage.dmp")
        output:
            os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        conda:
            "../envs/blast.yaml"
        log:
            "results/logs/finalData/BLAST.log"
        script:
            "../scripts/create_blast_taxonomy.py"

rule blast:
    input:
        "results/finalData/OTU_representatives.fasta",
        expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".nhr", ".nin", ".nog", ".nsd", ".nsi", ".nsq"] if config["blast"]["database"] == "SILVA" else "")
    output:
        temp("results/finalData/blast_taxonomy.tsv")
    threads: 150
    params:
        db_path=config["blast"]["db_path"],
        max_target_seqs=str(config["blast"]["max_target_seqs"]) if config["blast"]["database"] == "BLAST" else "1",
        ident=str(config["blast"]["ident"]),
        evalue=str(config["blast"]["evalue"]),
        out6=str(config["blast"]["out6"])
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -num_threads {threads} -query {input[0]} -db {params.db_path}"
        " -max_target_seqs {params.max_target_seqs}"
        " -perc_identity {params.ident} -evalue {params.evalue}"
        " -outfmt {params.out6} -out {output}"

if config["blast"]["database"] == "NCBI":

    rule ncbi_taxonomy:
        input:
            blast_result = "results/finalData/blast_taxonomy.tsv",
            lineage = os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        output:
            tax_lineage = temp("results/finalData/blast_taxonomic_lineage.tsv"),
            all_tax = "results/finalData/blast_taxonomy_all.tsv"
        params:
            max_target_seqs=config["blast"]["max_target_seqs"]
        conda:
            "../envs/blast.yaml"
        log:
            "results/logs/finalData/BLAST.log"
        script:
            "../scripts/ncbi_taxonomy.py"

elif config["blast"]["database"] == "SILVA":

    rule silva_taxonomy:
        input:
            blast_result = "results/finalData/blast_taxonomy.tsv",
            lineage = os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        output:
            temp("results/finalData/blast_taxonomic_lineage.tsv")
        conda:
            "../envs/blast.yaml"
        log:
            "results/logs/finalData/BLAST.log"
        script:
            "../scripts/silva_taxonomy.py"

rule merge_results:
    input:
        merged_swarm="results/finalData/swarm_table.csv",
        blast_result="results/finalData/blast_taxonomic_lineage.tsv"
    output:
        "results/finalData/filtered_blast_table.csv"
    conda:
        "../envs/merge_results.yaml"
    log:
        "results/logs/finalData/BLAST.log"
    script:
        "../scripts/merge_results.py"
