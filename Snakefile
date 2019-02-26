import pandas as pd

units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), 'fq2'])

if config['merge']['paired_End']:
    GROUP = [1,2]
else:
    GROUP = 1

rule all:
    input:
        'results/finalData/filtered_table.csv',
        'results/finalData/merged.swarms',
        'results/qc/multiqc_report.html' if config["general"]["multiqc"] else [],
        'results/finalData/figures/AmpliconDuo.RData',
        'results/finalData/filtered_blast_table.csv' if config["blast"]["blast"] else []

ruleorder: assembly > prinseq

include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
include: "rules/merging.smk"
include: "rules/clustering.smk"
include: "rules/blast.smk"