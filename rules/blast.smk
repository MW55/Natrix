if config["blast"]["database"] == "SILVA":

    rule make_silva_db:
        output:
            expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".nhr", ".nin", ".nog", ".nsd", ".nsi", ".nsq"])
        params:
            db_path=config["blast"]["db_path"]
        conda:
            "../envs/blast.yaml"
        shell:
            "dir_name=$(dirname {params[0]});"
            "wget -P $dir_name/ --progress=bar "
            "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/"
            "SILVA_132_SSURef_tax_silva.fasta.gz;"
            "gunzip $dir_name/SILVA_132_SSURef_tax_silva.fasta.gz;"
            "makeblastdb -in $dir_name/SILVA_132_SSURef_tax_silva.fasta "
            "-dbtype nucl -parse_seqids -out $dir_name/silva.db"

elif config["blast"]["database"] == "NCBI":

    rule make_ncbi_db:
        output:
            expand(config["blast"]["db_path"] + "*" + "{file_extension}", file_extension=[".nhr", ".nin", ".nog", ".nsd", ".nsi", ".nsq"])
        params:
            db_path=config["blast"]["db_path"]
        conda:
            "../envs/blast.yaml"
        shell:
            "dir_name=$(dirname {params[0]});"
            "wget -P $dir_name/ --progress=bar "
            "for i in {00..71}; "
            "do wget -P $dir_name/ --progress=bar "
            "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.${i}.tar.gz; "
            "tar xvzf $dir_name/nt.${i}.tar.gz -C $dir_name; "
            "done; touch 'nt'"

rule blast:
    input:
        "results/finalData/OTU_representatives.fasta",
        expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".nhr", ".nin", ".nog", ".nsd", ".nsi", ".nsq"] if config["blast"]["database"] == "SILVA" else "")
    output:
        temp("results/finalData/blast_taxonomy.tsv")
    threads: 150
    params:
        db_path=config["blast"]["db_path"],
        max_target_seqs=str(config["blast"]["max_target_seqs"]),
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

rule merge_results:
    input:
        merged="results/finalData/merged.swarms",
        final_table_path2="results/finalData/filtered_table.csv",
        blast_result="results/finalData/blast_taxonomy.tsv"
    output:
        "results/finalData/filtered_blast_table.csv"
    conda:
        "../envs/merge_results.yaml"
    log:
        "results/logs/finalData/BLAST.log"
    script:
        "../scripts/merge_results.py"
