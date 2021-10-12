import re

output_file=open(snakemake.output[0], "w")
input_file=open(snakemake.input[0], "r")

for line in input_file:
    arr = line.split(',')
    arr[-1] = re.sub(r"\([^()]*\)", "", arr[-1])
    str = ",".join(arr)
    output_file.write(str)
