import glob
import os
from snakemake.utils import validate

validate(config, "schema/config.schema.yaml")

SAMPLES = sorted(set([os.path.basename(x)[:-len('.fastq.gz')].split(config["sample"]["separator"])[config["sample"]["name_idx"]].replace('_','-') for x in glob.glob(config['general']['filename'] + '/*.fastq.gz')]))

if config["merge"]["filter_method"] == 'split_sample':
    UNITS = ['A','B']
else:
    UNITS = ['A']

if config["merge"]["paired_End"]:
    READS = [1,2]
else:
    READS = 1

rule all:
    input:
        "results/finalData/unfiltered_table.csv",
        "results/finalData/filtered_table.csv",
        "results/finalData/swarm_table.csv" if config["general"]["seq_rep"] == "OTU" else [],
        "results/qc/multiqc_report.html" if config["general"]["multiqc"] else [],
        "results/finalData/figures/AmpliconDuo.RData" if config["merge"]["ampliconduo"] and config["merge"]["filter_method"] == "split_sample" else [],
        "results/finalData/filtered_blast_table.csv" if config["blast"]["blast"] else [],
        "results/finalData/filtered_blast_table_complete.csv" if config["blast"]["blast"] else []


ruleorder: assembly > prinseq

include: "rules/preprocess_sample_names.smk"
include: "rules/demultiplexing.smk"
include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
include: "rules/merging.smk"
include: "rules/clustering.smk"
include: "rules/blast.smk"
