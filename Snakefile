import pandas as pd
from snakemake.utils import validate

validate(config, "schema/config.schema.yaml")

units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
name_ext = config["merge"]["name_ext"][:-1]

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), "fq2"])

if config["merge"]["paired_End"]:
    reads = [1,2]
else:
    reads = 1

ranks = ['Phylum','Class','Order','Family','Genus','Species']

rule all:
    input:
        "results/finalData/unfiltered_table.csv",
        "results/finalData/filtered_table.csv",
        "results/finalData/swarm_table.csv" if config["general"]["seq_rep"] == "OTU" else [],
        "results/qc/multiqc_report.html" if config["general"]["multiqc"] else [],
        "results/finalData/figures/AmpliconDuo.RData" if config["merge"]["ampliconduo"] and config["merge"]["filter_method"] == "split_sample" else [],
        "results/finalData/filtered_blast_table.csv" if config["blast"]["blast"] else [],
        "results/finalData/filtered_blast_table_complete.csv" if config["blast"]["blast"] else [],
        expand("results/finalData/taxonomy_splits/{rank}.csv", rank=ranks) if config["blast"]["blast"] and config["blast"]["split_filtered_blast_table"] else []


ruleorder: assembly > prinseq

include: "rules/demultiplexing.smk"
include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
include: "rules/merging.smk"
include: "rules/clustering.smk"
include: "rules/blast.smk"
