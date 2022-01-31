# this script is for generating fasta from mothur OTU Table file
import sys
import csv
f1 = open(snakemake.output[0], 'w')
with open(snakemake.input[0], 'r') as file:
        reader = csv.reader(file)
        next(reader, None)  # skip the headers
        for row in reader:
            f1.write('>' + row[0] + '\n' + row[1] + "\n")
