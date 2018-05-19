rule write_fasta:
    input:
        'results/finalData/filtered_table.csv'
    output:
        'results/finalData/filtered.fasta'
    run:
        import csv
        with open(input[0], 'r') as csv_in, open(output[0], 'w') as fasta_out:
                  filtered_table = csv.reader(csv_in)
                  next(filtered_table) # skip the header
                  for row in enumerate(filtered_table):
                      fasta_out.write('>{};size={};\n{}\n'.format(row[0],
                      sum([int(num) for num in row[1][1:]]), row[1][0]))

rule swarm:
    input:
        'results/finalData/filtered.fasta'
    output:
        'results/finalData/merged_representatives.fasta',
        'results/finalData/merged.swarms'
    threads: config['general']['cores']
    conda:
        '../envs/swarm.yaml'
    shell:
        'swarm -t {threads} -f -z -w {output[0]} < {input} > {output[1]}'
