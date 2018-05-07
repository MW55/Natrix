import dinopy
from glob import glob
from collections import Counter
import subprocess

#seqs = dinopy.FastaReader(str(snakemake.input))
#sorted_counter = Counter(seqs.reads()).most_common()
#sample_name = snakemake.wildcards.sample + '_' + snakemake.wildcards.unit

subprocess.call('cd-hit-est -i {} -o {}_cdhit.fasta -c {} -T {} -s {}'.format(snakemake.input, str(snakemake.input)[:-6], snakemake.params.id_percent, snakemake.params.cores, snakemake.params.length_cutoff), shell=True)

seq = dinopy.FastaReader(str(snakemake.input)[:-6] + '_cdhit.fasta')
clust = str(snakemake.input)[:-6] + '_cdhit.fasta.clstr'

def printcc(current, count):
        if current is not None and count > 0:
            return '{}; size = {};'.format(current.strip().replace('Cluster ',
                clust.split('/')[-2] + '_'), count)

def cluster_size(clust):
    clust_sizes = []
    with open(clust) as f:
        current = None
        count = 0
        for line in f:
            if line[0] == '>':
                clust_sizes.append(printcc(current, count))
                current = line
                count = 0
            else:
                count += 1
        clust_sizes.append(printcc(current, count))
    clust_sizes.pop(0)
    return clust_sizes

c_size = cluster_size(clust)
with dinopy.FastaWriter(str(snakemake.output), force_overwrite=True, line_width = 1000) as clust:
    for seq in enumerate(seqs.entries()):
        output.write_entry((seq[1].sequence, clust_sizes[seq[0]].encode()))
    output.close()

#with dinopy.FastaWriter(str(snakemake.output), line_width=1000) as clust:
#    for i in range(len(sorted_counter)):
#        name = '{}_{};size={};'.format(sample_name, i+1, sorted_counter[i][1])
#        clust.write_entry((sorted_counter[i][0],bytes(name, 'utf-8')))
#    clust.close()
