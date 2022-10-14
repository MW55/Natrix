import operator

import pandas as pd

skip_columns = ['sequences', 'qlen', 'length', 'pident', 'mismatch', 'qstart', 'qend', 'sstart', 'send', 'gaps',
                'evalue', 'seqid', 'Unnamed: 0']
blast_table = pd.read_csv(snakemake.input[0],  index_col='taxonomy', usecols=lambda x: x not in skip_columns)
if snakemake.params.database == 'NCBI':
    tax_keys = blast_table[['tax_key', "tax_ids"]]
    tax_keys = tax_keys.groupby(tax_keys.index).first()
    blast_table = blast_table.drop(['tax_key', 'tax_ids'], axis=1)
blast_table = blast_table.astype('int64')
#blast_table = blast_table.rename(columns=lambda x: x.split('-')[0])
blast_table = blast_table.T
blast_table = blast_table.groupby(by=blast_table.columns, axis=1).sum()

ranks = ['phylum','class','order','family','genus','species']

def rindex(lst, value):
    #actually returns index + 1
    return len(lst) - operator.indexOf(reversed(lst), value) if value in lst else -1

def split_taxonomy(x,rank):
    key = tax_keys.at[x,'tax_key']
    index = rindex(key.split(';'),rank)
    return x if index == -1 else ';'.join(x.split(';')[:index])

def split_ids(x,rank):
    key = tax_keys.at[x,'tax_key']
    index = rindex(key.split(';'),rank)
    return tax_keys.at[x,'tax_ids'].split(';')[-1] if index == -1 else tax_keys.at[x,'tax_ids'].split(';')[:index][-1]

for o in reversed(snakemake.output):
    rank = ranks.pop()
    if(snakemake.params.database == 'NCBI'):
        blast_table = blast_table.rename(columns=lambda x: split_taxonomy(x, rank))
        tax_keys = tax_keys.rename(index=lambda x: split_taxonomy(x, rank))
        tax_keys = tax_keys.groupby(tax_keys.index).first()

    blast_table = blast_table.groupby(by=blast_table.columns, axis=1).sum()
    if(snakemake.params.database == "NCBI"):
        id_row = pd.Series({col:split_ids(col, rank) for col in blast_table.columns}, name="tax_ids")
        blast_table = blast_table.append(id_row)
    blast_table.to_csv(o)
    if(snakemake.params.database == "NCBI"):
        blast_table = blast_table.drop(['tax_ids'], axis=0)

    if(snakemake.params.database == 'SILVA'):
        max_count = max([len(i.split(';')) for i in blast_table.columns])
        blast_table = blast_table.rename(columns=lambda x: x.rsplit(';', 1)[0] if len(x.split(';')) == max_count and len(x.split(';')) >= len(ranks) else x) 