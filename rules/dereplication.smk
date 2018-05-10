rule cdhit:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.fasta'
    output:
        'results/assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta'
    conda:
        '../envs/dereplication.yaml'
    params:
        id_percent = config['derep']['clustering'],
        cores = config['general']['cores'],
        length_cutoff = config['derep']['length_overlap']
    shell:'cd-hit-est -i {input} -o {output} -c {params.id_percent} -T' 
          '{params.cores} -s {params.length_cutoff} -M 2000 -sc 1 -d 0'

rule cluster_sorting:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta'
    output:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.fasta'
    conda:
        '../envs/dereplication.yaml'
    script:
        '../scripts/dereplication.py'
