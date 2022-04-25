import pandas as pd
import logging

logging.basicConfig(filename=str(snakemake.log),
            level=logging.DEBUG)

tax_df = pd.read_csv(snakemake.input[0], sep="\t\\|\t", index_col=0, header=None, engine='python',
                     names=['tax_id', 'lineage'])
nodes_df = pd.read_csv(snakemake.input[1], sep="\t\\|\t", index_col=0, header=None, engine='python', usecols=[0, 2],
                       names=['tax_id', 'rank'])
names_df = pd.read_csv(snakemake.input[2], sep="\t\\|\t", index_col=0, header=None, engine='python', usecols=[0, 1, 3],
                       names=['tax_id', 'name_txt', 'name class'])


def remove_endofline_tab(value):
    return value[:-2] if isinstance(value, str) else value


tax_df['lineage'] = tax_df['lineage'].apply(remove_endofline_tab)
names_df['name class'] = names_df['name class'].apply(remove_endofline_tab)

names_df = names_df[names_df['name class'] == 'scientific name']
names_df = names_df.drop(columns=['name class'])
names_df = names_df.reindex(names_df.index.values.sort())
names_df = names_df.join(nodes_df)


def ids_to_names_and_ranks(ids):
    ids = list(map(int, ids.split(" ")))
    names = [names_df.at[id, 'name_txt'] for id in ids]
    ranks = [names_df.at[id, 'rank'] for id in ids]
    return ";".join(names), ";".join(ranks)


def get_tax_lineage(row):
    if not isinstance(row['lineage'], str):
        names, ranks = ids_to_names_and_ranks(str(row.name))
        logging.info("Taxid {} is not present in NCBI lineage file.".format(row.name))  # no lineage information?
    else:
        names, ranks = ids_to_names_and_ranks(row['lineage'] + str(row.name))
    return names, ranks

tax_df[['tax_lineage', 'tax_key']] = tax_df.apply(get_tax_lineage, axis=1, result_type="expand")
tax_df = tax_df.drop(['lineage'], axis=1)

tax_df.to_hdf(snakemake.output[0], key='tax_df', mode='w')
