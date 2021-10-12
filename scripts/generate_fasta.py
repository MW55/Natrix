# this script is for generating fasta from mothur OTU Table file
import sys
import csv
f1 = open(snakemake.output[0], 'w')
with open(snakemake.input[0], 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            f1.write('>' + row[0] + '\n' + row[1] + "\n")

# this block is to remove headers(first two lines)  from fasta file
with open(snakemake.output[0], 'r+') as fp:
    # read an store all lines into list
    lines = fp.readlines()
    # move file pointer to the beginning of a file
    fp.seek(0)
    # truncate the file
    fp.truncate()

    # start writing lines
    # iterate line and line number
    for number, line in enumerate(lines):
        # delete line number 5 and 8
        # note: list index start from 0
        if number not in [0, 1]:
            fp.write(line)
