rule blast:
    input:
        'results/finalData/merged_representatives.fasta'
    output:
        'results/finalData/blast_taxonomy.tsv'
    threads:
        config['general']['cores']
    params:
        db_path = config['blast']['db_path'],
        max_target_seqs = str(config['blast']['max_target_seqs']),
        ident = str(config['blast']['ident']),
        evalue = str(config['blast']['evalue']),
        out6 = str(config['blast']['out6'])
    conda:
        '../envs/blast.yaml'
    shell:
        'blastn -num_threads {threads} -query {input} -db {params.db_path}'
        ' -max_target_seqs {params.max_target_seqs} -perc_identity {params.ident}'
        ' -evalue {params.evalue} -outfmt {params.out6} -out {output}'
