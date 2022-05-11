import pandas as pd
import collections as col

filtered_table = pd.read_csv(snakemake.input["final_table_path2"],
                index_col="seqid")
columns = filtered_table.columns
filtered_table = filtered_table.to_dict(orient="index")

representatives = dict()

with open(snakemake.input["consensus"], "r") as rep_in:
    centroid = ""
    sequence = ""
    for line in rep_in:
        if line.startswith(">"):
            representatives[centroid] = sequence
            header = line.split(";")
            if(snakemake.params["consensus"] == "consensus"):
                centroid = header[0].split("=")[1]+";"+header[1]+";"
            else:
                centroid = header[0][1:]+";"+header[1]+";"
            sequence = ""
        else:
            sequence += line.rstrip()
    representatives[centroid] = sequence

swarm_dict = dict()
filtered = dict()
with open(snakemake.input["merged"], "r") as swarms:
    for row in swarms:
        splitted_row = [s for s in row.split()]
        if (splitted_row[0] == "C"):
            cluster_id = splitted_row[8]
            cluster_num = splitted_row[2]
            if cluster_id not in swarm_dict:
                swarm_dict[cluster_id] = dict.fromkeys(columns, 0)
            swarm_dict[cluster_id]["seqid"] = "N{}_{}".format(cluster_id.split(";")[0],cluster_num)
            swarm_dict[cluster_id]["sequences"] = representatives[cluster_id]
            for sample in filter(lambda i: i not in ["sequences", "seqid"], columns):
                swarm_dict[cluster_id][sample] += filtered_table[">" + cluster_id][sample]
        elif splitted_row[0] == "H":
            cluster_id = splitted_row[9]
            query_id = splitted_row[8]
            if cluster_id not in swarm_dict:
                swarm_dict[cluster_id] = dict.fromkeys(columns, 0)
            for sample in filter(lambda i: i not in ["sequences", "seqid"], columns):
                swarm_dict[cluster_id][sample] += filtered_table[">" + query_id][sample]

representatives_filtered = ""
for k, v in swarm_dict.items():
    if sum([value for sample, value in swarm_dict[k].items() if not sample in ["sequences", "seqid"]]) > snakemake.params["min_sequences"]:
        filtered[k] = v
        representatives_filtered += ">"+v["seqid"]+"\n"
        representatives_filtered += "\n".join([v["sequences"][i: i+60] for i in range(0, len(v["sequences"]), 60)])+"\n"

df = pd.DataFrame.from_dict(filtered, orient="index").set_index("seqid")
df.to_csv(snakemake.output["swarm_table"])

df = pd.DataFrame.from_dict(swarm_dict, orient="index").set_index("seqid")
df.to_csv(snakemake.output["all_out"])

with open(snakemake.output["consensus_filtered"], "w") as repr_out:
    repr_out.write(representatives_filtered)
