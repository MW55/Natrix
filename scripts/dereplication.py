import dinopy
from glob import glob
from collections import Counter

seqs = dinopy.FastaReader(str(snakemake.input))
sorted_counter = Counter(seqs.reads()).most_common()
sample_name = snakemake.wildcards.sample + '_' + snakemake.wildcards.unit
with dinopy.FastaWriter(str(snakemake.output), line_width=1000) as clust:
    for i in range(len(sorted_counter)):
        name = '{}_{};size={};'.format(sample_name, i+1, sorted_counter[i][1])
        clust.write_entry((sorted_counter[i][0],bytes(name, 'utf-8')))
    clust.close()
