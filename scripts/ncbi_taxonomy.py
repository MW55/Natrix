import pandas as pd
import logging
import math

columns = ["seqid","qlen","length","pident","mismatch","qstart","qend",
            "sstart","send","gaps","evalue","taxonomy", "dbseqid"]

logging.basicConfig(filename=str(snakemake.log),
            level=logging.DEBUG)

blast_df = pd.read_csv(snakemake.input["blast_result"], sep="\t",names=columns, index_col="seqid")

taxids = list(set(blast_df["taxonomy"]))

tax_df = pd.read_hdf(snakemake.input["lineage"], 'tax_df')
tax_df = tax_df.reindex(taxids)

def split2(char, string):
    if(isinstance(string, str)):
        if char in string:
            return string.split(char)
        else:
            return [string]
    else:
        return ["NaN"]


def create_filter_set():
    res_set = set()
    for f in snakemake.params.drop_tax_classes.split(","):
        f = f.strip()
        if not f:
            pass
        elif f.isdigit():
            res_set.add(int(f))
        else:
            # f = f.replace("*", ".*")
            tmp = tax_df[tax_df.tax_lineage.str.contains(f, na=False)]
            res_set.update(tmp.index)
    return res_set

querynames = list(set(blast_df.index.values))

header = True
f_set = create_filter_set()
for i in querynames:
    blast_res = blast_df.loc[i,:] # blast results
    if blast_res.ndim == 1: # case only 1 blast hit
        blast_res = blast_res.to_frame().T
    blast_res = blast_res[~blast_res['taxonomy'].isin(f_set)] # remove blast results with blacklisted taxonomy
    tax_res = tax_df.loc[blast_res['taxonomy'],:] # taxonomy for result
    tax_res = tax_res.set_index(blast_res.index)
    blast_res = pd.concat([blast_res, tax_res], axis=1)
    if blast_res.empty:
        continue
    blast_res["tax_level"] = blast_res.apply(lambda row: math.log(min(len(split2("; ", row['tax_lineage'])),6)), axis=1) # max length 6, log
    blast_res["len_val"] = blast_res.apply(lambda row: math.log(row['qlen']+1), axis=1)
    blast_res["Mscore"] = blast_res.apply(lambda row: row['pident']*row['tax_level']*row['len_val'], axis=1)
    blast_res = blast_res.sort_values(by=['Mscore'], ascending=False) # sort by score
    mode = 'w' if header else 'a'
    blast_res.to_csv(snakemake.output["all_tax"], mode=mode, header=header, sep="\t", index_label="seqid")
    best_tax = blast_res[:1][["qlen","length","pident","mismatch","qstart","qend", "sstart","send","gaps","evalue","tax_lineage", "tax_key"]]
    best_tax = best_tax.rename(index=str, columns={"tax_lineage": "taxonomy"})
    best_tax.to_csv(snakemake.output["tax_lineage"], mode=mode, header=header, sep="\t", index_label="seqid")
    header = False
