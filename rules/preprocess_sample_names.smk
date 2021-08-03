rule preprocess_sample_names:
    output:
        temp(expand("input/{sample}_{unit}_R{read}.fastq.gz", sample=SAMPLES,unit=UNITS,read=READS)),
        temp("primertable.csv.tmp")
    params:
        sample_config=config['sample'],
        filename_base=config['general']['filename'],
        primertable=config['general']['primertable'],
        paired=config['merge']['paired_End'],
        filter_method=config['merge']['paired_End']
    log:
        "results/logs/preprocess_sample_names.log"
    run:
        import glob
        import os
        import logging
        import pandas as pd

        logging.basicConfig(filename=str(log),level=logging.DEBUG)

        samples = glob.glob(params.filename_base+ '/*.fastq.gz')
        logging.info(f"""Processing {len(samples)} Samples with following configuration:
separator: {params.sample_config['separator']}
name_idx: {params.sample_config['name_idx']}
unit_idx: {params.sample_config['unit_idx']}
read_idx: {params.sample_config['read_idx']}
unit_A_identifier: {params.sample_config['unit_A_identifier']}
read_forward_identifier: {params.sample_config['read_forward_identifier']}""")
        suffix = '.fastq.gz'

        if not os.path.exists('input/'):
            os.mkdir('input')

        primertable = pd.read_csv(params.primertable, index_col="Probe")

        def map_sample(sample):
            splitted = os.path.basename(sample)[:-len(suffix)].split(params.sample_config['separator'])
            mapped = splitted[params.sample_config['name_idx']].replace('_','-')
            mapped += '_'
            mapped += 'A' if params.filter_method=='not_split' or params.sample_config['unit_A_identifier'] == splitted[params.sample_config['unit_idx']] else 'B'
            mapped += '_'
            mapped += 'R1' if not params.paired or params.sample_config['read_forward_identifier'] == splitted[params.sample_config['read_idx']] else 'R2'
            mapped += suffix
            return mapped

        mapped_samples = {s: map_sample(s) for s in samples}

        if(len(mapped_samples.values()) != len(set(mapped_samples.values()))):
            logging.warning("Some files are mapped to the same name!")
        sample_len = len(os.path.basename(max(mapped_samples.keys(),key=lambda f:len(os.path.basename(f)))))
        logging.info("Files are mapped as follows:\n"+"\n".join(f"{os.path.basename(old_file).ljust(sample_len)} => {new_file}" for old_file, new_file in mapped_samples.items()))

        for old_file,new_file in mapped_samples.items():
            os.symlink(os.path.abspath(old_file),'input/' + new_file)

        def map_primertable(sample):
            splitted = sample.split(params.sample_config['separator'])
            mapped = splitted[params.sample_config['name_idx']].replace('_','-')
            mapped += '_'
            mapped += 'A' if params.filter_method == 'not_split' or params.sample_config['unit_A_identifier'] == splitted[params.sample_config['unit_idx']] else 'B'
            return mapped

        primertable.rename(index=lambda x: map_primertable(x),inplace=True) #change primertable to multiindex?
        primertable.to_csv("primertable.csv.tmp")