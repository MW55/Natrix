import dinopy

fasta = snakemake.input[0]
clstr = snakemake.input[1]
fasta_not_clstr = snakemake.input[2]

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
def not_clstr_dict(clstr_list, seq_dict_mc):
    not_clstr = dict()
    for id in clstr_list:
        s = seq_dict_mc[id.encode()]
        if s in not_clstr.keys():
            not_clstr[s]["names"].append(id)
            not_clstr[s]["count"] += 1
        else:
            not_clstr[s] = {"count" : 1}
            not_clstr[s]["names"] = list()
            not_clstr[s]["names"].append(id)
    return not_clstr

def count_clstr(clstr_list, seq_dict_mc):
    not_clstr = not_clstr_dict(clstr_list, seq_dict_mc)
    count = 0
    sum = 0
    max = None
    for s in not_clstr.keys():
        if not_clstr[s]["count"] > count:
            max = s
            count = not_clstr[s]["count"]
        sum += not_clstr[s]["count"]
    return [sum, not_clstr[max]["names"][0]]

# Version that uses the most common sequence as representative,
# because of the larger overhead it is much slower than taking
# the longest sequence as representative, but could lead to
# a better representation of the genetic composition of the
# sample.
def get_most_common_rep(clstr, fasta_not_clstr):
    seq_dict_mc = sequence_dict(fasta_not_clstr)
    with open(clstr) as f:
        clstr_list = []
        cluster = -1
        clust_sizes = []
        ids = []
        for line in f:
            if line[0] == ">":
                if cluster != -1:
                    if snakemake.params.length_cutoff == 0:
                        clust_sizes.append("{};size={};".format(clstr.split("/")[-2] + "_" + str(cluster), str(len(clstr_list))))
                        ids.append(clstr_list[0])
                    else:
                        hd = count_clstr(clstr_list, seq_dict_mc)
                        clust_sizes.append("{};size={};".format(clstr.split("/")[-2] + "_" + str(cluster), str(hd[0])))
                        ids.append(hd[1])
                clstr_list = []
                cluster += 1
            else:
                clstr_list.append(line[line.find(">")+1:line.find("...")])
        # for last cluster
        if snakemake.params.length_cutoff == 0:
            clust_sizes.append("{};size={};".format(clstr.split("/")[-2] + "_" + str(cluster), str(len(clstr_list))))
            ids.append(clstr_list[0])
        else:
            hd = count_clstr(clstr_list, seq_dict_mc)
            clust_sizes.append("{};size={};".format(clstr.split("/")[-2] + "_" + str(cluster), str(hd[0])))
            ids.append(hd[1])
    return list(zip(clust_sizes, ids))

# Writes the representatives in a fasta file with the new headers, uses
# the old header as key to get the matching sequences from the dict.
def writer(c_size, seq_dict):
    with dinopy.FastaWriter(str(snakemake.output), force_overwrite=True,
            line_width=1000) as clust:
        clust.write_entries([(seq_dict[line[1].encode()],
            line[0].encode()) for line in c_size])
        clust.close()

if str(snakemake.params.repr) == "most_common":
    seq_dict_mc = sequence_dict(fasta_not_clstr)
    c_size_mc = get_most_common_rep(clstr, fasta_not_clstr)
    writer(c_size_mc, seq_dict_mc)
else:
    seq_dict_l = sequence_dict(fasta)
    c_size_l = get_longest_rep(clstr)
    writer(c_size_l, seq_dict_l)
