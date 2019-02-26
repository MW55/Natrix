rule write_fasta:
    input:
        'results/finalData/filtered_table_temp.csv'
    output:
        'results/finalData/filtered.fasta',
        'results/finalData/filtered_table.csv'
    run:
        import csv
        with open(input[0], 'r') as csv_in, open(
            output[0], 'w') as fasta_out, open(output[1], 'w') as csv_out:
            filtered_table = csv.reader(csv_in)
            filtered_table_seqid = csv.writer(csv_out)
            for row in enumerate(filtered_table):
                if row[0] != 0:
                    abu_sum = sum([int(num) for num in row[1][1:]])
                    filtered_table_seqid.writerow(['>{};size={};'.format(
                        row[0], abu_sum)] + row[1])
                    fasta_out.write('>{};size={};\n{}\n'.format(row[0],
                        abu_sum, row[1][0]))
                else:
                    filtered_table_seqid.writerow(['seqid'] + row[1])

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
