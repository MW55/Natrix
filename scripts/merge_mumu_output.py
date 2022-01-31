import pandas as pd
import re

output_all=open(snakemake.output[0], "w")
output_table=open(snakemake.output[1], "w")
output_meta=open(snakemake.output[2], "w")
mumu_file=open(snakemake.input[0], "r")
otu_file=open(snakemake.input[1], "r")

data=pd.read_csv(otu_file, index_col=0)
mumu=pd.read_csv(mumu_file, sep='\t', index_col=0)

data_inter=data.loc[mumu.index]

out=pd.concat([data_inter["sequences"], mumu, data_inter["taxonomy"]], axis=1)

out.to_csv(output_all)
mumu.to_csv(output_table)

out2 = data_inter[["sequences", "taxonomy"]]

# substitute bootstrap values if present and append to dataframe
taxonomy2 = data_inter["taxonomy"].tolist()
if re.match("(\d*)", taxonomy2[0]):
    taxonomy2 = list(map(lambda x: re.sub(r"\([^()]*\)", "", x), taxonomy2))
    pd.options.mode.chained_assignment = None
    out2.loc[:,'taxonomy2'] = taxonomy2

out2.to_csv(output_meta)
