import pandaseq as pd
import yaml

with open(str(snakemake.input), 'r') as f_:
    seq_dict = yaml.load(f_)

#### no split sample: filter : cutoff filter

### split sample(A+B): use the func that compares the length of a list to a set
