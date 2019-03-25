import dinopy
from glob import glob

fasta = snakemake.input[0]
clstr = snakemake.input[1]
fasta_not_clstr = snakemake.input[2]

print(str(snakemake.params))

def sequence_dict(fasta_file):
    return {entry.name:entry.sequence for entry
                in dinopy.FastaReader(str(fasta_file)).entries()}


# Returns a list of tuples containing the old header of the
# representative seq of each cluster and the new header as required
# by uchime (sorted by cluster size).
# This version uses the longest sequence as representative of a cluster,
# just as CDHIT does.
def get_longest_rep(clstr):
    clust_sizes = []
    ids = []
    with open(clstr) as f:
        current = None
        count = 0
        for line in f:
            if line[0] == ">":
                if current is not None and count > 0:
                    clust_sizes.append("{};size={};".format(
                        current[1:].strip().replace("Cluster ",
                            clstr.split("/")[-2] + "_"), count))
                current = line
                count = 0
            elif line[-2] == "*":
                ids.append(line[line.find(">")+1:line.find("...")])
                count += 1
            else:
                count += 1
    return list(zip(clust_sizes, ids))

# Helper functions to get the abundance of each sequence per
# sequence header, needed to use the most common sequence
# of a cluster as representative.
def not_clstr_dict(not_clstr_fasta):
    print("not_clstr_dict")
    not_clstr = dict()
    for entry in dinopy.FastaReader(not_clstr_fasta).entries():
        s = entry.sequence.decode()
        if s in not_clstr.keys():
            not_clstr[s]["names"].add(entry.name.decode())
            not_clstr[s]["count"] += 1
        else:
            not_clstr[s] = {"count" : 1}
            not_clstr[s]["names"] = set()
            not_clstr[s]["names"].add(entry.name.decode())
    return not_clstr

def header_dict(not_clstr_fasta):
    print("not_clstr_fasta")
    not_clstr = not_clstr_dict(not_clstr_fasta)
    return {n:not_clstr[s]["count"] for s 
            in not_clstr.keys() for n 
            in not_clstr[s]["names"]}

def dictfilt(header_dict, seq_set):
    return {h : header_dict[h] for h in header_dict if h in seq_set}

# Version that uses the most common sequence as representative,
# because of the larger overhead it is much slower than taking
# the longest sequence as representative, but could lead to
# a better representation of the genetic composition of the
# sample.
def get_most_common_rep(clstr, fasta_unclstr):
    print("get_most_common_seq")
    clust_sizes = []
    ids = []
    seqs = dict()
    hd = header_dict(fasta_unclstr)
    with open(clstr) as f:
        current = None
        count = 0
        for line in f:
            if line[0] == ">":
                if current is not None and count > 0:
                    clust_sizes.append("{};size={};".format(
                        current[1:].strip().replace("Cluster ",
                            clstr.split("/")[-2] + "_"), count))
                    ids.append(max(dictfilt(hd, seqs[current])))
                current = line
                count = 0
                seqs[current] = set()
            elif line[-2] == "*":
                seqs[current].add(line[line.find(">")+1:line.find("...")])
                count += 1    
            else:
                seqs[current].add(line[line.find(">")+1:line.find("...")])
                count += 1
    return list(zip(clust_sizes, ids))

# Writes the representatives in a fasta file with the new headers, uses
# the old header as key to get the matching sequences from the dict.
def writer(c_size, seq_dict):
    print("writer")
    with dinopy.FastaWriter(str(snakemake.output), force_overwrite=True,
            line_width=1000) as clust:
        clust.write_entries([(seq_dict[line[1].encode()],
            line[0].encode()) for line in c_size])
        clust.close()

if str(snakemake.params) == "most_common":
    seq_dict_mc = sequence_dict(fasta_not_clstr)
    c_size_mc = get_most_common_rep(clstr, fasta_not_clstr)
    writer(c_size_mc, seq_dict_mc)
else:
    seq_dict_l = sequence_dict(fasta)
    c_size_l = get_longest_rep(clstr)
    writer(c_size_l, seq_dict_l)