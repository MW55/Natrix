def get_fastq(wildcards, read_pair='fq1'):
    return units.loc[(wildcards.sample, wildcards.unit), [read_pair]].dropna()
rule unzip:
    input:
         data_folder + '/{sample}_{unit}_R{read}.fastq.gz'
    output:
        data_folder + '/{sample}_{unit}_R{read}.fastq'
    shell: 'gunzip -c {input} > {output}'

#das define_primer script muss noch für single end sets angepasst werden, das die auch functionieren!
rule define_primer:
    input:
        primer_table = config['general']['filename'] + '.csv'
    output:
        'p_table_with_primer.csv'
    script:
        '../scripts/define_primer.R'

rule prinseq:
    input:
        r1 = lambda wildcards: data_folder + '/' + get_fastq(wildcards, 'fq1'), r2 = lambda wildcards: data_folder + '/' + get_fastq(wildcards, 'fq2')
    output:
        'results/assembly/{sample}_{unit}/{sample}_{unit}_1.fastq','results/assembly/{sample}_{unit}/{sample}_{unit}_2.fastq'
    params:
        config['qc']['mq']
    run:
        output_edit = str(output[0])[:-8]
        output_bad = str(output[0])[:-8] + '_bad'
        shell('bin/prinseq-lite/prinseq-lite.pl -verbose -fastq {input.r1} -fastq2 {input.r2} -ns_max_n 0 -min_qual_mean {params} -out_good {output_edit} -out_bad {output_bad} 2>&1')

if config["merge"]["paired_End"]:
    rule assembly:
        input:
            r1 = 'results/assembly/{sample}_{unit}/{sample}_{unit}_1.fastq',
            r2 = 'results/assembly/{sample}_{unit}/{sample}_{unit}_2.fastq',
            primer_t = 'p_table_with_primer.csv'
        output:
             'results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq'
        params:
            threshold = config["qc"]["threshold"],
            minoverlap = config["qc"]["minoverlap"],
            minlen = config["qc"]["minlen"],
            maxlen = config["qc"]["maxlen"],
            minqual = config["qc"]["minqual"],
            prim_rm = config["qc"]["all_primer"]
        log:
            'logs/read_assembly.log'
        run:
            primer_table = pd.read_csv(input.primer_t)
            r1_primer = primer_table['f_primer'].loc[primer_table['Probe_FWD'] == wildcards.sample + '_' + wildcards.unit].values[0]
            r2_primer = primer_table['r_primer'].loc[primer_table['Probe_FWD'] == wildcards.sample + '_' + wildcards.unit].values[0]
            if params.prim_rm:
                shell('pandaseq -f {input.r1} -r {input.r2} -B -a -F -g {log} -w {output} -N -t  {params.threshold}, -o {params.minoverlap} -l {params.minlen} -L {params.maxlen} -C min_phred:{params.minqual}')
            else:
                shell('pandaseq -f {input.r1} -r {input.r2} -B -a -F -g {log} -w {output} -N -p "{r1_primer}" -q "{r2_primer}" -t  {params.threshold}, -o {params.minoverlap} -l {params.minlen} -L {params.maxlen} -C min_phred:{params.minqual}')
else:
    rule digDeeper:
        input:
            file_ = data_folder + '{sample}_{unit}_R1.fastq',
            file_name = '{sample}_{unit}'
        output:
        #got to find a dataset&primertable for single end and adjust the digDeeper script
            'results/digDeeper_done'
        run:
             r1_primer = primer_table['r_primer'].loc[primer_table['Probe_FWD'] == wildcards.sample + '_' + wildcards.unit].values[0],
             primer = '../script/apply.DigDeeper1.03.R'

rule copy_to_fasta:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq'
    output:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.fasta'
    shell:
        'cat {input} | paste - - - - | cut -f 1,2 | sed "s/^@/>/g" | tr "\t" "\n" > {output}'
