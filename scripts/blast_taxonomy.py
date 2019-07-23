import pandas as pd
import logging
import math

columns = ["seqid","qlen","length","pident","mismatch","qstart","qend",
            "sstart","send","gaps","evalue","taxonomy"]

logging.basicConfig(filename=str(snakemake.log),
            level=logging.DEBUG)

blast_df = pd.read_csv(snakemake.input["blast_result"], sep="\t",names=columns, index_col="seqid")

taxids = list(set(blast_df["taxonomy"]))

tax_df = pd.read_hdf(snakemake.input["lineage"], 'tax_df')
tax_df = tax_df.loc[taxids]

def split2(char, string):
    if(isinstance(string, str)):
        if char in string:
            return string.split(char)
        else:
            return [string]
    else:
        return ["NaN"]

querynames=list(set(blast_df.index.values))

header = True
for i in querynames:
    blast_res = blast_df.loc[i,:] # blast results
    if blast_res.ndim == 1: # case only 1 blast hit
        blast_res = blast_res.to_frame().T
    tax_res = tax_df.loc[blast_res['taxonomy'],:] # taxonomy for result
    tax_res = tax_res.set_index(blast_res.index)
    blast_res = pd.concat([blast_res, tax_res], axis=1)
    blast_res["tax_level"] = blast_res.apply(lambda row: math.log(min(len(split2("; ", row['tax_lineage'])),6)), axis=1) # max length 6, log
    blast_res["len_val"] = blast_res.apply(lambda row: math.log(row['qlen']+1), axis=1)
    blast_res["Mscore"] = blast_res.apply(lambda row: row['pident']*row['tax_level']*row['len_val'], axis=1)
    blast_res = blast_res.sort_values(by=['Mscore'], ascending=False) # sort by score
    mode = 'w' if header else 'a'
    blast_res.to_csv(snakemake.output["all_tax"], mode=mode, header=header, sep="\t")
    best_tax = blast_res[:1][["qlen","length","pident","mismatch","qstart","qend", "sstart","send","gaps","evalue","tax_lineage"]]
    best_tax = best_tax.rename(index=str, columns={"tax_lineage": "taxonomy"})
    best_tax.to_csv(snakemake.output["tax_lineage"], mode=mode, header=header, sep="\t")
    header = False
