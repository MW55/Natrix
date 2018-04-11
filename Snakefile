import pandas as pd

configfile: "rev_test.yaml"

# Create .tsv files containing the samples and units in an ordered fashion
file_list = os.listdir(config['general']['filename'])
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
        "logs/qc_done"
       # "results/qc/multiqc_report.html"

include: "rules/quality_control.smk"
