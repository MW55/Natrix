import pandas as pd
import glob

configfile: "test_data.yaml"

# Create .tsv files containing the samples and units in an ordered fashion
data_folder = config['general']['filename']
file_list = [file.split('/')[-1] for file in glob.glob(data_folder + "/*.gz")]
file_list.sort()
df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
    index =range(int(len(file_list)/2)))
if config['merge']['paired_End']:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
        index =range(int(len(file_list)/2)))
    i, j = (0, 0)
    while i < len(file_list)/2:
        df.loc[i] = file_list[j].split('_')[0], file_list[j].split('_')[1], \
            file_list[j][:-3], file_list[j+1][:-3]
        j += 2
        i += 1
    reads = [1,2]
else:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
        index = range(int(len(file_list))))
    i = 0
    while i < len(file_list):
        df.loc[i] = file_list[i].split('_')[0], file_list[i].split('_')[1], \
            file_list[i][:-3], np.nan
        i += 1
    reads = 1
df.to_csv('units.tsv', sep='\t')
pd.DataFrame(df['sample']).to_csv('samples.tsv', sep='\t')
samples = pd.read_table(config["general"]["samples"], index_col="sample")
units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

rule all:
    input: 'results/qc/multiqc_report.html', expand([data_folder + '/{unit.sample}_{unit.unit}_R{read}.fastq','results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_{read}.fastq','results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_assembled.fastq','results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}.fasta', 'logs/qc_done', 'results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}.clustered100.fasta'], unit = units.reset_index().itertuples(), read=reads), 'logs/qc_done'


include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
