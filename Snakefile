import pandas as pd
from snakemake.utils import validate

#validate(config, "schema/config.schema.yaml")

units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
name_ext = config["merge"]["name_ext"][:-1]

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), "fq2"])

def get_fastq(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        return expand("demultiplexed/{sample}_{unit}_{group}.fastq",
                        group=[1,2], **wildcards)
    return "demultiplexed/{sample}_{unit}_1.fastq".format(**wildcards)

if config["merge"]["paired_End"]:
    reads = [1,2]
else:
    reads = 1

"""
rule all:
    input:
        "results/finalData/unfiltered_table.csv" if not config["shotgun"]["enabled"] else [],
        "results/finalData/filtered_table.csv" if not config["shotgun"]["enabled"] else [],
        "results/finalData/swarm_table.csv" if config["general"]["seq_rep"] == "OTU" and not config["shotgun"]["enabled"] else [],
        "results/qc/multiqc_report.html" if config["general"]["multiqc"] else [],
        "results/finalData/figures/AmpliconDuo.RData" if config["merge"]["ampliconduo"] and config["merge"]["filter_method"] == "split_sample" else [],
        "results/finalData/filtered_blast_table.csv" if config["blast"]["blast"] else [],
        "results/finalData/filtered_blast_table_complete.csv" if config["blast"]["blast"] else [],
        #"filtered/{unit.sample}_clean_R1.fastq.gz" if config["shotgun"]["enabled"] else [],
        #"filtered/{unit.sample}_clean_R2.fastq.gz" if config["shotgun"]["enabled"] else [],
        "results/taxonomy/taxonomy_abundance_matrix.tsv" if config["shotgun"]["enabled"] else [],
        "results/functional/functional_abundance_matrix.tsv" if config["shotgun"]["enabled"] else [],
        "results/diversity/alpha_diversity.tsv" if config["shotgun"]["enabled"] else [],
        "results/diversity/beta_diversity_pcoa.tsv" if config["shotgun"]["enabled"] else [],
        "results/visualization/taxonomy_barplot.png" if config["shotgun"]["enabled"] else [],
        "results/visualization/functional_heatmap.png" if config["shotgun"]["enabled"] else [],
        "results/visualization/beta_diversity_pcoa.png" if config["shotgun"]["enabled"] else []
"""

rule all:
    input:
        "results/rgi/merged_rgi.txt"
        #"results/rgi/{unit.sample}_{unit.unit}_rgi.txt"


ruleorder: assembly > prinseq

include: "rules/demultiplexing.smk"
include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
include: "rules/merging.smk"
include: "rules/clustering.smk"
include: "rules/blast.smk"
include: "rules/shotgun_metagenome/host_removal.smk"
include: "rules/shotgun_metagenome/read_based_processing.smk"
include: "rules/shotgun_metagenome/functional_annotation.smk"
include: "rules/shotgun_metagenome/aggregate_results.smk"
include: "rules/shotgun_metagenome/analysis.smk"
include: "rules/antimicrobial_resistance_prediction/resistance_prediction.smk"