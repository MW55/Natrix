rule dereplication:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.fasta'
    output:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.fasta'
    conda:
        '../envs/dereplication.yaml'
    params:
        id_percent = config['derep']['clustering'],
        cores = config['general']['cores'],
        length_cutoff = config['derep']['length_overlap']
    script:
        '../scripts/dereplication.py'
