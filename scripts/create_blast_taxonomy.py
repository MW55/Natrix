import pandas as pd
import logging

logging.basicConfig(filename=str(snakemake.log),
            level=logging.DEBUG)

tax_df = pd.read_csv(snakemake.input[0], sep="\t\\|\t", index_col=0, header=None, engine='python')

def get_tax_lineage(row):
    if(isinstance(row[2], str)):
        result = row[2][:-2]
    else:
        result = row[2]
        logging.info("Taxid {} is not present in NCBI lineage file.".format(row.name))
    return result

tax_df["tax_lineage"] = tax_df.apply(lambda row: get_tax_lineage(row), axis=1).map(str) + tax_df[1]
tax_df = tax_df.drop([1,2], axis=1)

tax_df.to_hdf(snakemake.output[0], key='tax_df', mode='w')
