rule unfiltered_table:
    input:
        expand('results/finalData/{unit.sample}_{unit.unit}.clustered100.nonchimera.fasta', unit=units.reset_index().itertuples())
    output:
        'results/finalData/unfiltered_table.csv',
        temp('results/finalData/unfiltered_dict.yaml')
    conda:
        '../envs/unfiltered_table.yaml'
    script:
        '../scripts/unfiltered_table.py'

rule filtering:
    input:
        'results/finalData/unfiltered_dict.yaml'
    output:
        'results/finalData/filtered_table.csv',
        'results/finalData/filtered_out_table.csv',
        temp('results/finalData/filtered_dict.yaml'),
        temp('results/finalData/filtered_out_dict.yaml')
    params:
        filter_method = config['merge']['filter_method'],
        cutoff = config['merge']['cutoff']
    conda:
        '../envs/filtering.yaml'
    script:
        '../scripts/filtering.py'
