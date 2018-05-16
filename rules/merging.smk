rule unfiltered_table:
    input:
        expand('results/finalData/{unit.sample}_{unit.unit}.clustered100.nonchimera.fasta', unit=units.reset_index().itertuples())
    output:
        'results/finalData/unfiltered_table.csv',
        'results/finalData/unfiltered_dict.json' #temp()
    conda:
        '../envs/unfiltered_table.yaml'
    script:
        '../scripts/unfiltered_table.py'

rule filtering:
    input:
        'results/finalData/unfiltered_dict.json'
    output:
        'results/finalData/filtered_table.csv',
        'results/finalData/filtered_out_table.csv',
        temp('results/finalData/filtered_dict.json'),
        temp('results/finalData/filtered_out_dict.json')
    params:
        filter_method = config['merge']['filter_method'],
        cutoff = config['merge']['cutoff']
    conda:
        '../envs/filtering.yaml'
    script:
        '../scripts/filtering.py'
