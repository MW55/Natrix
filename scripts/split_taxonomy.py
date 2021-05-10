import operator

import pandas as pd

skip_columns = ['sequences', 'qlen', 'length', 'pident', 'mismatch', 'qstart', 'qend', 'sstart', 'send', 'gaps',
                'evalue', 'seqid', 'Unnamed: 0']
blast_table = pd.read_csv(snakemake.input[0], index_col='taxonomy', usecols=lambda x: x not in skip_columns)

tax_keys = blast_table['tax_key'] if snakemake.params.database == 'NCBI' else None
tax_keys = tax_keys.T
blast_table = blast_table.astype('int64')
blast_table = blast_table.rename(columns=lambda x: x.split('-')[0])
blast_table = blast_table.T


ranks = ['phylum','class','order','family','genus','species']

def rindex(lst, value):
    #actually returns index + 1
    return len(lst) - operator.indexOf(reversed(lst), value) if value in lst else -1

def split_taxonomy(x,rank):
    key = tax_keys.at['tax_key',x]
    index = rindex(key.split(';'),rank)
    return x if index == -1 else x[:index]

for o in reversed(snakemake.output):
    rank = ranks.pop()
    if(snakemake.params.database == 'NCBI'):
        blast_table = blast_table.rename(columns=lambda x: split_taxonomy(x, rank))
        tax_keys = tax_keys.rename(columns=lambda x: split_taxonomy(x, rank))
        tax_keys = tax_keys.groupby(by=tax_keys.columns, axis=1).first()

    blast_table = blast_table.groupby(by=blast_table.columns, axis=1).sum()
    blast_table.to_csv(o)

    if(snakemake.params.database == 'SILVA'):
        max_count = max([len(i.split(';')) for i in blast_table.columns])
        blast_table = blast_table.rename(columns=lambda x: x.rsplit(';', 1)[0] if len(x.split(';')) == max_count and len(x.split(';')) >= len(ranks) else x)