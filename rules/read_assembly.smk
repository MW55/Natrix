
rule define_primer:
    input:
        primer_table = config['general']['filename'] + '.csv'
    output:
        'p_table_with_primer.csv'
    script:
        '../scripts/define_primer.R'

rule prinseq:
    input:
        prinseq = 'bin/prinseq-lite/prinseq-lite.pl',
       # fastq1 = data_folder + '/{sample}_{unit}_R1.fastq',
       # fastq2 = data_folder + '/{sample}_{unit}_R2.fastq'
 
        fastq1 = expand(data_folder + '/{unit.sample}_{unit.unit}_R1.fastq', unit = units.reset_index().itertuples(), group=groups),
        fastq2 = expand(data_folder + '/{unit.sample}_{unit.unit}_R2.fastq', unit = units.reset_index().itertuples(), group=groups)
    output:
        #'results/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_R{group}.fastq'
        fq_good = expand('results/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_R{group}.fastq',
        unit = units.reset_index().itertuples(), group=groups),
        fq_bad = expand('results/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_R{group}_bad.fastq',
        unit = units.reset_index().itertuples(), group=groups)
    params:
        config['qc']['mq']
    shell:
        '{input.prinseq} -verbose -fastq {input.fastq1} -fastq2 {input.fastq2}'
        ' -ns_max_n 0 -min_qual_mean {params} -out_good {output.fq_good}'
        ' -out_bad {output.fq_bad}_bad 2>&1'
        
#rule assembly:
#    input:
#        data_folder = config['general']['filename'],
#        primer_table = config['general']['filename'] + '.csv',
#    output:
#        "logs/assembly_done"
#        #dir = '/results/assembly'
#    log:
#        "logs/read_assembly.log",
#	"logs/pandaseq.log"
#    script:
#        "../scripts/read_assembly.R"
