import dinopy
from glob import glob

seq_dict = {entry.name:entry.sequence for entry
            in dinopy.FastaReader(str(snakemake.input[0])).entries()}

# Returns a list of tuples containing the old header of the
# representative seq of each cluster and the new header as required
# by uchime (sorted by cluster size).
clust_sizes = []
ids = []
with open(snakemake.input[1]) as f:
    current = None
    count = 0
    for line in f:
        if line[0] == ">":
            if current is not None and count > 0:
                clust_sizes.append("{};size={};".format(
                    current[1:].strip().replace("Cluster ",
                        snakemake.input[1].split("/")[-2] + "_"), count))
            current = line
            count = 0
        elif line[-2] == "*":
            ids.append(line[line.find(">")+1:line.find("...")])
            count += 1
        else:
            count += 1
c_size = list(zip(clust_sizes, ids))

# Writes the representatives in a fasta file with the new headers, uses
# the old header as key to get the matching sequences from the dict.
with dinopy.FastaWriter(str(snakemake.output), force_overwrite=True,
        line_width=1000) as clust:
    clust.write_entries([(seq_dict[line[1].encode()],
        line[0].encode()) for line in c_size])
    clust.close()
