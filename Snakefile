import pandas as pd
import os
from snakemake.utils import validate

validate(config, "schema/config.schema.yaml")

units = pd.read_table(os.path.join(config["general"]["output_dir"],config["general"]["units"]), index_col=["sample", "unit"],
    dtype=str)
#print(units)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
name_ext = config["merge"]["name_ext"][:-1]
#print(units.index)

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), "fq2"])

if config["merge"]["paired_End"]:
    reads = [1,2]
else:
    reads = 1

# rules
rule all:
    input:
        os.path.join(config["general"]["output_dir"],"qc/multiqc_report.html") if config["general"]["multiqc"] else [],
        os.path.join(config["general"]["output_dir"],"filtering/unfiltered_table.csv"),
        os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
        os.path.join(config["general"]["output_dir"],"filtering/figures/AmpliconDuo.RData") if config["merge"]["ampliconduo"] and config["merge"]["filter_method"] == "split_sample" else [],
        os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" and config['dataset']['nanopore'] == 'FALSE' else [],
        #mothur
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table.csv"), database=config['classify']['database']) if config['classify']['mothur'] else [],
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.csv"), database=config['classify']['database']) if config["general"]["seq_rep"] == "OTU" and  config['classify']['mothur'] and config['postcluster']['mumu'] else [],
        # blast
        expand(os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/OTU_table.csv"), database=config['blast']['database'].lower()) if config["blast"]["blast"] and config['dataset']['nanopore'] == 'FALSE' else [],
        expand(os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/OTU_table_mumu.csv"), database=config['blast']['database'].lower()) if config["general"]["seq_rep"] == "OTU" and  config["blast"]["blast"] and config['postcluster']['mumu'] else [],
        # nanopore

if not config['dataset']['nanopore']:
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
#include: "rules/denoise.smk"
include: "rules/pychop.smk"
include: "rules/read_correction.smk"
