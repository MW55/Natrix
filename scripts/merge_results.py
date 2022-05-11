import pandas as pd
import collections as col
import logging
import numpy as np

blast = pd.read_csv(snakemake.input['blast_result'], sep="\t")

if snakemake.params['seq_rep'] == 'ASV':
    merged = pd.read_csv(snakemake.input['merged'], sep=",") #merged = pd.read_csv(snakemake.input['merged'], index_col=0, sep=",")
else:
    merged = pd.read_csv(snakemake.input['merged'], index_col=0, sep=",")


def new_index(table):
    new_index = table["seqid"].tolist()[0:]
    if snakemake.params["platform"] == "Illumina":
        splitted = [i.lstrip(">").replace("size=", "").split(";")[0:2] for i in new_index]
        new_index = ["N{}_{}".format(i[0], i[1]) for i in splitted]
    else:
        new_index = [i.lstrip(">") for i in new_index]
    df = pd.DataFrame(new_index)
    table["seqid"] = df[0]
    table = table.set_index('seqid')
    return table

blast = new_index(blast)
if snakemake.params['seq_rep'] == 'ASV':
    merged = new_index(merged)

result = merged.join(blast, how='outer')

logging.basicConfig(filename=str(snakemake.log), level=logging.DEBUG)

nas = pd.isna(result["taxonomy"])
nbh_counter = sum(nas)
bh_counter = sum(~nas)

no_blast_hit = list(result.loc[nas, :].index)

logging.info("{} sequences could be assigned taxonomic information using BLAST,\
            {} could not be assigned taxonomic information using BLAST".format(bh_counter, nbh_counter))
logging.info("No fitting BLAST hit found for seq_id:")
logging.info(no_blast_hit)

cols = set(result.columns.tolist())
cols_include = ['sequences', 'qlen', 'length', 'pident', 'mismatch', 'qstart', 'qend', 'sstart', 'send', 'gaps', 'evalue', 'taxonomy']
if snakemake.params['database'] == 'NCBI':
    cols_include.append('tax_key')
    cols_include.append('tax_ids')
cols_all = cols_include + list(cols - set(cols_include))

result = result[cols_all]
result.to_csv(snakemake.output["complete"], index_label="seqid")
# filter out OTUs without taxonomic annotation
result = result.loc[~nas, :]
result.to_csv(snakemake.output["filtered"], index_label="seqid")
