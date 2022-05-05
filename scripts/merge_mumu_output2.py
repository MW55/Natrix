import pandas as pd
import re

mumu_file=open(snakemake.input[0], "r")
otu_file=open(snakemake.input[1], "r")

output_all=open(snakemake.output[0], "w")
output_table=open(snakemake.output[1], "w")
output_meta=open(snakemake.output[2], "w")

data=pd.read_csv(otu_file, index_col=0)
mumu=pd.read_csv(mumu_file, sep='\t', index_col=0)

data_inter=data.loc[mumu.index]

out=pd.concat([data_inter[["sequences", "qlen", "length", "pident", "mismatch", "qstart", "qend", "sstart", "send", "gaps", "evalue"]], mumu, data_inter["taxonomy"]], axis=1)

out.to_csv(output_all)
mumu.to_csv(output_table)

out2 = data_inter[["sequences", "qlen", "length", "pident", "mismatch", "qstart", "qend", "sstart", "send", "gaps", "evalue", "taxonomy"]]
out2.to_csv(output_meta)
