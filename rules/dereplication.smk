rule dereplication:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.fasta'
    output:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.fasta'
    conda:
        '../envs/dereplication.yaml'
    script:
        '../scripts/dereplication.py'
