
if config["classify"]["database"] == "PR2":
        rule download_pr2:
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
elif config["classify"]["database"] == "UNITE":
        rule download_unite:
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

elif config["classify"]["database"] == "SILVA":
        rule download_silva:
            output:
                expand(["database/silva_db.{silva_db_version}.fasta", "database/silva_db.{silva_db_version}.tax"], silva_db_version=config["database_version"]["silva"])
            params:
                db_version=config["database_version"]["silva"],
                scripts=config["scripts"]["silva"] #to access scripts in scripts directory
            shell:
                """
                wget -P ./ --progress=bar -O silva_{params.db_version}.fasta.gz https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_{params.db_version}_SSURef_tax_silva.fasta.gz;
                gunzip -c silva_{params.db_version}.fasta.gz > database/silva_db.{params.db_version}.fasta;
                wget -P ./ --progress=bar -O silva_{params.db_version}.tax.gz https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_{params.db_version}.txt.gz
                gunzip -c silva_{params.db_version}.tax.gz > database/silva_db.{params.db_version}.tax;
                bash {params.scripts}/silva_db_edit.sh database/silva_db.{params.db_version}.tax > database/silva.tax.tmp; #to convert silva taxonomy file to mothur format
                rm database/silva_db.{params.db_version}.tax;
                mv database/silva.tax.tmp database/silva_db.{params.db_version}.tax;
                sed 1d -i database/silva_db.{params.db_version}.tax #to remove header in taxonomy file
                """

