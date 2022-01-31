
if config["classify"]["database"] == "pr2":
        rule download_pr2:
            input:
                expand(os.path.join(config["general"]["output_dir"],"finalData/representatives.fasta"))
            output:
                expand(["database/pr2db.{pr2_db_version}.fasta", "database/pr2db.{pr2_db_version}.tax"], pr2_db_version=config["database_version"]["pr2"])
            params:
                db_version=config["database_version"]["pr2"]
            shell:
                """
                wget -P ./ --progress=bar https://github.com/pr2database/pr2database/releases/download/v{params.db_version}/pr2_version_{params.db_version}_SSU_mothur.fasta.gz;
                wget -P ./ --progress=bar https://github.com/pr2database/pr2database/releases/download/v{params.db_version}/pr2_version_{params.db_version}_SSU_mothur.tax.gz;
                gunzip -c ./pr2_version_{params.db_version}_SSU_mothur.fasta.gz > database/pr2db.{params.db_version}.fasta;
                gunzip -c ./pr2_version_{params.db_version}_SSU_mothur.tax.gz > database/pr2db.{params.db_version}.tax;
                rm ./pr2_version_{params.db_version}_SSU_mothur.tax.gz ./pr2_version_{params.db_version}_SSU_mothur.fasta.gz;
                """
elif config["classify"]["database"] == "unite":
        rule download_unite:
            input:
                expand(os.path.join(config["general"]["output_dir"],"finalData/representatives.fasta"))
            output:
                expand("database/unite_v8.3.fasta"), expand("database/unite_v8.3.tax")
            shell:
                """
                wget -P ./ --progress=bar  -O UNITE_public_mothur_10.05.2021.tgz --progress=bar https://files.plutof.ut.ee/public/orig/38/6A/386A46113D04602A78FB02497D9B0E1A8FE2145B23C2A6314A62B419F0D08E73.tgz;
                tar -xvf UNITE_public_mothur_10.05.2021.tgz;
                mv UNITE_public_mothur_10.05.2021/UNITE_public_mothur_10.05.2021_taxonomy.txt database/unite_v8.3.tax;
                mv UNITE_public_mothur_10.05.2021/UNITE_public_mothur_10.05.2021.fasta {output[0]};
                rm -rf UNITE_public_mothur_10.05.2021;
                """

elif config["classify"]["database"] == "silva":

        rule download_silva:
            input:
                expand(os.path.join(config["general"]["output_dir"],"finalData/representatives.fasta"))
            output:
                expand(["database/silva_db.{silva_db_version}.fasta", "database/silva_db.{silva_db_version}.tax.temp"], silva_db_version=config["database_version"]["silva"])
            params:
                db_version=config["database_version"]["silva"]
            shell:
                """
                wget -P ./ --progress=bar -O database/silva_{params.db_version}.fasta.gz https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_{params.db_version}_SSURef_tax_silva.fasta.gz;
                gunzip -c database/silva_{params.db_version}.fasta.gz > database/silva_db.{params.db_version}.fasta;
                wget -P ./ --progress=bar -O database/silva_{params.db_version}.tax.gz https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_{params.db_version}.txt.gz
                gunzip -c database/silva_{params.db_version}.tax.gz > database/silva_db.{params.db_version}.tax.temp;
                """

        rule edit_silva:
           input:
              expand("database/silva_db.{silva_db_version}.tax.temp", silva_db_version=config["database_version"]["silva"])
           output:
              expand("database/silva_db.{silva_db_version}.tax", silva_db_version=config["database_version"]["silva"])   
           conda:
              "../envs/blast.yaml"
           script:
              "../scripts/edit_silva_mothur.py"
