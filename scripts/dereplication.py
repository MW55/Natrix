from glob import glob

fasta = snakemake.input[0]
clstr = snakemake.input[1]
fasta_not_clstr = snakemake.input[2]

# FASTA parser
def read_fasta(file):
    sequences = {}
    with open(file, 'r') as f:
        header = None
        seq = []
        for line in f:
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(seq)
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if header:
            sequences[header] = ''.join(seq)
    return sequences

# FASTA writer
def write_fasta(file, records):
    with open(file, 'w') as f:
        for header, sequence in records.items():
            f.write(f">{header}\n{sequence}\n")

def sequence_dict(fasta_file):
    return read_fasta(fasta_file)

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

def not_clstr_dict(clstr_list, seq_dict_mc):
    not_clstr = {}
    for id in clstr_list:
        s = seq_dict_mc[id]
        if s in not_clstr:
            not_clstr[s]["names"].append(id)
            not_clstr[s]["count"] += 1
        else:
            not_clstr[s] = {"count": 1, "names": [id]}
    return not_clstr

def count_clstr(clstr_list, seq_dict_mc):
    not_clstr = not_clstr_dict(clstr_list, seq_dict_mc)
    count = 0
    max_seq = None
    for s, data in not_clstr.items():
        if data["count"] > count:
            max_seq = s
            count = data["count"]
    return [sum(data["count"] for data in not_clstr.values()), not_clstr[max_seq]["names"][0]]

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
        if snakemake.params.length_cutoff == 0:
            clust_sizes.append("{};size={};".format(clstr.split("/")[-2] + "_" + str(cluster), str(len(clstr_list))))
            ids.append(clstr_list[0])
        else:
            hd = count_clstr(clstr_list, seq_dict_mc)
            clust_sizes.append("{};size={};".format(clstr.split("/")[-2] + "_" + str(cluster), str(hd[0])))
            ids.append(hd[1])
    return list(zip(clust_sizes, ids))

def writer(c_size, seq_dict, output_file):
    with open(output_file, 'w') as f:
        for new_header, old_header in c_size:
            sequence = seq_dict.get(old_header, "")
            f.write(f">{new_header}\n{sequence}\n")

if str(snakemake.params.repr) == "most_common":
    seq_dict_mc = sequence_dict(fasta_not_clstr)
    c_size_mc = get_most_common_rep(clstr, fasta_not_clstr)
    writer(c_size_mc, seq_dict_mc, snakemake.output[0])
else:
    seq_dict_l = sequence_dict(fasta)
    c_size_l = get_longest_rep(clstr)
    writer(c_size_l, seq_dict_l, snakemake.output[0])
