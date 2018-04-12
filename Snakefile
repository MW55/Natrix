import pandas as pd

configfile: "test_data.yaml"

# Create .tsv files containing the samples and units in an ordered fashion
data_folder = config['general']['filename']
file_list = os.listdir(data_folder)
file_list.sort()
df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
	index =range(int(len(file_list)/2)))

if config['merge']['paired_End']:
	df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
		index =range(int(len(file_list)/2)))
	i, j = (0, 0)
	while i < len(file_list)/2:
		df.loc[i] = file_list[j].split('_')[0], file_list[j].split('_')[1], \
			file_list[j], file_list[j+1]
		j += 2
		i += 1
 	groups = [1,2]	
else:
	df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
		index = range(int(len(file_list))))
	i = 0
	while i < len(file_list):
		df.loc[i] = file_list[i].split('_')[0], file_list[i].split('_')[1], \
			file_list[i], np.nan
		i += 1
	groups = 1
df.to_csv('units.tsv', sep='\t')
pd.DataFrame(df['sample']).to_csv('samples.tsv', sep='\t')

samples = pd.read_table(config["general"]["samples"], index_col="sample")
units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],	
	dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

rule all:
    input:
         expand('results/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_R{group}.fastq',
         unit = units.reset_index().itertuples(), group = groups)


       # "results/qc/multiqc_report.html"
        # expand('results/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_R{group}.fastq',
       # unit = units.reset_index().itertuples(), group = groups)

       # "p_table_with_primer.csv"
       # expand(data_folder +
       # '/{unit.sample}_{unit.unit}_R{group}.fastq',
       # unit = units.reset_index().itertuples(), group = groups),
         # "logs/assembly_done"
       # "logs/qc_done"

rule unzip:
    input:
        data_folder + '/{sample}_{unit}_R{group}.fastq.gz'
    output:
        data_folder + '/{sample}_{unit}_R{group}.fastq'
    shell:
        'for i in {input}; do (gzip -kd $i) done'

include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
