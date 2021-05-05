rule preprocess_sample_names:
    output:
        temp(expand("input/{sample}_{unit}_R{read}.fastq.gz", sample=SAMPLES,unit=UNITS,read=READS)),
        temp("primertable.csv.tmp")
    params:
        sample_config=config['sample'],
        filename_base=config['general']['filename'],
        primertable=config['general']['primertable'],
        paired=config['merge']['paired_end'],
        filter_method=config['merge']['paired_end']
    run:
        import glob
        import os
        import pandas as pd

        samples = glob.glob(params.filename_base+ '/*.fastq.gz')
        suffix = '.fastq.gz'

        if not os.path.exists('input/'):
            os.mkdir('input')

        primertable = pd.read_csv(params.primertable, index_col="Probe")

        for file in samples:
            splitted = os.path.basename(file)[:-len(suffix)].split(params.sample_config['separator'])
            mapped = splitted[params.sample_config['name_idx']].replace('_','-')
            mapped += '_'
            mapped += 'A' if params.filter_method=='not_split' or params.sample_config['unit_A_identifier'] == splitted[params.sample_config['unit_idx']] else 'B'

            primertable.rename(index=lambda x: mapped if x in os.path.basename(file) else x, inplace=True) #probably replace this

            mapped += '_'
            mapped += 'R1' if not params.paired or params.sample_config['read_forward_identifier'] == splitted[params.sample_config['read_idx']] else 'R2'

            mapped += suffix

            os.symlink(os.path.abspath(file),'input/'+mapped)

        primertable.to_csv("primertable.csv.tmp")