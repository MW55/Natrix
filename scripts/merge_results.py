import pandas as pd
import collections as col
import logging
import numpy as np

blast = pd.read_csv(snakemake.input['blast_result'], sep="\t")

swarm = pd.read_csv(snakemake.input['merged_swarm'], index_col=0, sep=",")

new_index = blast["seqid"].tolist()[0:]
splitted = [i.replace("size=", "").split(";")[0:2] for i in new_index]
new_index = ["N{}_{}".format(i[0], i[1]) for i in splitted]

df = pd.DataFrame(new_index)
blast["seqid"] = df[0]

blast = blast.set_index('seqid')
result = swarm.join(blast, how='outer')

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
cols_all = cols_include + list(cols - set(cols_include))

result = result.loc[~nas, cols_all]
result.to_csv(snakemake.output[0], index_label="seqid")
