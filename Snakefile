import pandas as pd
import os
from snakemake.utils import validate

validate(config, "schema/config.schema.yaml")

units = pd.read_table(os.path.join(config["general"]["output_dir"],config["general"]["units"]), index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
name_ext = config["merge"]["name_ext"][:-1]

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), "fq2"])

if config["merge"]["paired_End"]:
    reads = [1,2]
else:
    reads = 1

rule all:
    input:
        os.path.join(config["general"]["output_dir"],"finalData/unfiltered_table.csv"),
        os.path.join(config["general"]["output_dir"],"finalData/filtered_table.csv"),
        os.path.join(config["general"]["output_dir"],"finalData/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" else [],
        os.path.join(config["general"]["output_dir"],"qc/multiqc_report.html") if config["general"]["multiqc"] else [],
        os.path.join(config["general"]["output_dir"],"finalData/figures/AmpliconDuo.RData") if config["merge"]["ampliconduo"] and config["merge"]["filter_method"] == "split_sample" else [],
        #mothur
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/mothur_out.taxonomy"), database=config['classify']['database']) if config['classify']['mothur'] else [],
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/swarm_mothur.csv"), database=config['classify']['database']) if config['classify']['mothur'] else [],
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.csv"), database=config['classify']['database']) if config['mumu']['mumu'] else [],
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/FINAL_OUTPUT_OTU.csv"), database=config['classify']['database']) if config['mumu']['mumu'] else [],
        "database/silva_db.138.1.fasta" if config["classify"]["database"] == "silva" else [], 
        "database/silva_db.138.1.tax.temp" if config["classify"]["database"] == "silva" else [],
#       expand("database/silva_db.{silva_db_version}.tax", silva_db_version=config["database_version"]["silva"]) if config["classify"]["database"] == "silva" else [],
        "database/pr2db.4.14.0.fasta" if config["classify"]["database"] == "pr2" else [],
        "database/unite_v8.3.fasta" if config["classify"]["database"] == "unite" else [], 
        # blast
#       expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/filtered_blast_table.csv"), database=config['blast']['output']) if config["blast"]["blast"] else [],
#       expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/filtered_blast_table_complete.csv"), database=config['blast']['output']) if config["blast"]["blast"] else []
        os.path.join(config["general"]["output_dir"],"finalData/blast_silva/filtered_blast_table.csv") if config["blast"]["database"] == "SILVA"else [],
        os.path.join(config["general"]["output_dir"],"finalData/blast_silva/filtered_blast_table_complete.csv") if config["blast"]["database"] == "SILVA"else [],
        os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/filtered_blast_table_complete.csv") if config["blast"]["database"] = "NCBI" else [],
        os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/filtered_blast_table.csv") if config["blast"]["database"] =="NCBI" else []

ruleorder: assembly > prinseq

include: "rules/demultiplexing.smk"
include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
include: "rules/merging.smk"
include: "rules/clustering.smk"
include: "rules/blast.smk" 
include: "rules/pr2_unite_silva.smk" 
include: "rules/classify.smk" 
include: "rules/mumu.smk"
