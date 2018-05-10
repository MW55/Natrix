import dinopy
from glob import glob
from collections import OrderedDict

seq_dict = {entry.name:entry.sequence for entry in dinopy.FastaReader(snakemake.input).entries()}
clust = str(snakemake.input)[:-6] + '_cdhit.fasta.clstr'

def cluster_size(clust):
    clust_sizes = []
    ids = []
    with open(clust) as f:
        current = None
        count = 0
        for line in f:
            if line[0] == '>':
                if current is not None and count > 0:
                    clust_sizes.append('{}; size = {};'.format(current[1:].strip().replace('Cluster ', clust.split('/')[-2] + '_'), count))
                current = line
                count = 0
            elif line[-2] == '*':
                ids.append(line[line.find('>')+1:line.find('...')])
                count += 1
            else:
                count += 1
        clust_sizes.append(printcc(current, count))
    return list(zip(clust_sizes, ids))

c_size = cluster_size(clust)
with dinopy.FastaWriter(str(snakemake.output), force_overwrite=True, line_width = 1000) as clust:
    clust.write_entries([(sequence_dict[line[1].encode()], line[0].encode()) for line in c_size])
    clust.close()

#with dinopy.FastaWriter(str(snakemake.output), line_width=1000) as clust:
#    for i in range(len(sorted_counter)):
#        name = '{}_{};size={};'.format(sample_name, i+1, sorted_counter[i][1])
#        clust.write_entry((sorted_counter[i][0],bytes(name, 'utf-8')))
#    clust.close()
