import pandas as pd
import collections as col

filtered_dict = pd.read_csv(snakemake.input["final_table_path2"],
                index_col="seqid").to_dict(orient="index")

with open(snakemake.input["merged"], "r") as swarms:
    seq_names = []
    for row in swarms:
        otu_seq_list = [s for s in row.split()]
        seq_names.append(otu_seq_list)

swarm_dict = col.defaultdict()
for j in range(len(seq_names)):
    swarm_dict[j] = {}
    swarm_dict[j]["sequences"] = filtered_dict[">" + seq_names[j][0]]["sequences"]
    for sample in filter(lambda i: i not in ["sequences"], filtered_dict[">"
                    + seq_names[j][0]].keys()):
        swarm_dict[j][sample] = sum([filtered_dict[">"
            + key][sample] for key in seq_names[j]])
    swarm_dict[j]["seqid"] = "N{}_{}".format(seq_names[j][0].split(";")[0],
        sum([value for key, value in swarm_dict[j].items() if not key in {"seqid", "sequences"}]))

df = pd.DataFrame.from_dict(swarm_dict, orient="index").set_index("seqid")
df.to_csv(snakemake.output[0])