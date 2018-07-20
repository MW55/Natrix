import pandas as pd
import collections as col

columns = ['seqid','qlen','length','pident','mismatch','qstart','qend',
            'sstart','send','gaps','evalue','taxonomy']
blast_dict = pd.read_csv(snakemake.input['blast_result'],
                sep='\t',names = columns, index_col='seqid').to_dict(
                        orient = 'index')
blast_dict = dict((k.split(';')[0], v) for k, v in blast_dict.items())
filtered_dict = pd.read_csv(snakemake.input['final_table_path2'],
                index_col='seqid').to_dict(orient = 'index')

with open(snakemake.input['blast_result'], 'r') as swarms:
    seq_names = []
    for row in swarms:
        otu_seq_list = [s for s in row.split()]
        seq_names.append(otu_seq_list)

swarm_dict = col.defaultdict()
for j in range(len(seq_names)):
    swarm_dict[j] = {}
    swarm_dict[j]['sequences'] = filtered_dict['>' + seq_names[j][0]]['sequences']
    for sample in filter(lambda i: i!='sequences', filtered_dict['>'
                    + seq_names[j][0]].keys()):
        swarm_dict[j][sample] = sum([filtered_dict['>'
            + key][sample] for key in seq_names[j]])
    swarm_dict[j]['seqid'] = 'N{}_{}'.format(seq_names[j][0].split(';')[0],
            sum([value for key, value in swarm_dict[j].items() if key != 'sequences']))
    try:
        swarm_dict[j].update(blast_dict[seq_names[j][0].split(';')[0]])
    except KeyError:
        pass #put in some logging command, telling the user which and how many sequences
            # were without blasthit

df = pd.DataFrame.from_dict(swarm_dict, orient='index').set_index('seqid')
df.to_csv(snakemake.output[0])
