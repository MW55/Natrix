import pandas as pd
import logging
import dinopy

far = dinopy.FastaReader(snakemake.input[0])
header = [i.name.decode().split(" ", 1) for i in far.entries()]
df = pd.DataFrame(header, columns=["id", "taxonomy"])
df = df.set_index(keys="id", drop=True)

df.to_hdf(snakemake.output[0], key='df', mode='w')
