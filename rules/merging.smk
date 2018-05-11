rule unfiltered_table:
    input:
        expand('results/finalData/{sample}_{unit}.clustered100.nonchimera.fasta')
    output:
        'results/finalData/unfiltered_table.csv'
    env:
        '../envs/merging.yaml'
    script:
        '../scripts/unfiltered_table.py'
