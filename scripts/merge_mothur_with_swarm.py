import sys
import pandas as pd
import re
taxonomy = pd.read_csv(sys.argv[2], sep='\t', names=['seqid', 'taxonomy']) #read mothur taxonomy file
#print(taxonomy.head())
swarm_count = pd.read_csv(sys.argv[1]) #read swarm output file
print (type(taxonomy))
tax=taxonomy['taxonomy']
ID=taxonomy['seqid'].replace(r"^", "N", regex=True).replace(r"\;size\=", "_", regex=True).replace(r"\;", "", regex=True) #change seqid as per swarm seqid
concat_col=pd.concat([ID, tax], axis=1) #concatenate the id and taxonomy column
final_file=pd.merge(swarm_count, concat_col, left_on='seqid', right_on='seqid', how='left') #merge mothur taxonomy to swarm abundance file
print(final_file.head())
final_file.to_csv(sys.argv[3], index=False) #output file
