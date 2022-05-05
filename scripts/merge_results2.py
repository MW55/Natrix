import pandas as pd
import re

otu_file=open(snakemake.input[1], "r")
taxonomy=open(snakemake.input[0], "r")

output_all=open(snakemake.output[0], "w")
output_table=open(snakemake.output[1], "w")
output_meta=open(snakemake.output[2], "w")

taxonomy = pd.read_csv(taxonomy, sep='\t', names=['seqid', 'taxonomy']) # read mothur taxonomy file
#print(taxonomy.head())
otu_count = pd.read_csv(otu_file) #read otu table

tax=taxonomy['taxonomy']
ID=taxonomy['seqid'].replace(r"^", "N", regex=True).replace(r"\;size\=", "_", regex=True).replace(r"\;", "", regex=True) # change seqid as per otu seqid
concat_col=pd.concat([ID, tax], axis=1) #concatenate the id and taxonomy column

if(otu_count['seqid'][1][0] == ">"):
	otu_count['seqid']=otu_count['seqid'].replace(r"^>", "N", regex=True).replace(r"\;size\=", "_", regex=True).replace(r"\;", "", regex=True) # change ID if ASV table

final_file=pd.merge(otu_count, concat_col, left_on='seqid', right_on='seqid', how='left') # merge mothur taxonomy to otu abundance file

final_file.to_csv(output_all, index=False) #output file

out2 = final_file[["seqid", "sequences", "taxonomy"]]

# substitute bootstrap values if present and append to dataframe
taxonomy2 = out2["taxonomy"].tolist()
if re.match("(\d*)", taxonomy2[0]):
    taxonomy2 = list(map(lambda x: x if pd.isna(x) else re.sub(r"\([^()]*\)", "", x), taxonomy2))
    pd.options.mode.chained_assignment = None
    out2.loc[:,'taxonomy2'] = taxonomy2

out2.to_csv(output_meta, index=False)

# get only counts
out3=otu_count.drop(["sequences"], axis=1)
out3.to_csv(output_table, index=False)
