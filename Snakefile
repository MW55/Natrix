import pandas as pd

configfile: "test_data.yaml"

samples = pd.read_table(config["general"]["samples"], index_col="sample")
units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), 'fq2'])

if config['merge']['paired_End']:
    reads = [1,2]
else:
    reads = 1

rule all:
    input:
       # expand(['demultiplexed/{unit.sample}_{unit.unit}_{read}.fastq',
       # 'results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_{read}.fastq',
       # 'results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_assembled.fastq',
       # 'results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}.fasta',
       # 'results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}.clustered100.fasta',
       # 'results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}.clustered100.fasta',
       # 'results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}.clustered100.uchime.txt'],
       # unit = units.reset_index().itertuples(), read=reads),
       # 'results/qc/multiqc_report.html',
       # 'results/finalData/unfiltered_table.csv',
       # 'results/finalData/filtered_table.csv',
        'results/finalData/figures/AmpliconDuo.RData',
        'results/finalData/merged.swarms',
        'results/finalData/blast_taxonomy.tsv'

ruleorder: assembly > prinseq

# 'logs/qc_done'

include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
include: "rules/merging.smk"
include: "rules/clustering.smk"
include: "rules/blast.smk"
