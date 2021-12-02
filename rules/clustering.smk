# dada runs per sample after preprocessing and before chimera removal and AmpliconDuo
rule DADA2:
    input:
        fwd = expand(
            os.path.join(config["general"]["output_dir"],"results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_1_cut.fastq"), unit=units.reset_index().itertuples()),
        rev = expand(
            os.path.join(config["general"]["output_dir"],"results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_2_cut.fastq"), unit=units.reset_index().itertuples()) if config["merge"]["paired_End"] == True else []
    output:
        expand(os.path.join(config["general"]["output_dir"],"results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_dada.fasta"), unit=units.reset_index().itertuples())
    params:
        paired_end=config["merge"]["paired_End"],
        minoverlap=config["qc"]["minoverlap"],
        splitsamples=config["merge"]["filter_method"]
    conda:
        "../envs/dada2.yaml"
    log:
        os.path.join(config["general"]["output_dir"],"logs/dada2.log")
    script:
        "../scripts/dada2.R"

# SWARM clustering runs on all samples after AmpliconDuo
rule write_fasta:
    input:
        os.path.join(config["general"]["output_dir"],"results/finalData/filtered_table_temp.csv")
    output:
        os.path.join(config["general"]["output_dir"],"results/finalData/filtered.fasta"),
        os.path.join(config["general"]["output_dir"],"results/finalData/filtered_table.csv")
    run:
        import csv
        with open(input[0], "r") as csv_in, open(
            output[0], "w") as fasta_out, open(output[1], "w") as csv_out:
            filtered_table = csv.reader(csv_in)
            filtered_table_seqid = csv.writer(csv_out)
            for row in enumerate(filtered_table):
                if row[0] != 0:
                    abu_sum = sum([int(num) for num in row[1][1:]])
                    filtered_table_seqid.writerow([">{};size={};".format(
                        row[0], abu_sum)] + row[1])
                    fasta_out.write(">{};size={};\n{}\n".format(row[0],
                        abu_sum, row[1][0]))
                else:
                    filtered_table_seqid.writerow(["seqid"] + row[1])

rule swarm:
    input:
        os.path.join(config["general"]["output_dir"],"results/finalData/filtered.fasta")
    output:
        os.path.join(config["general"]["output_dir"],"results/finalData/representatives.fasta"),
        temp(os.path.join(config["general"]["output_dir"],"results/finalData/merged.swarms"))
    threads: config["general"]["cores"]
    conda:
        "../envs/swarm.yaml"
    shell:
        "swarm -t {threads} -f -z -w {output[0]} < {input} > {output[1]}"

rule swarm_results:
    input:
        merged=os.path.join(config["general"]["output_dir"],"results/finalData/merged.swarms"),
        final_table_path2=os.path.join(config["general"]["output_dir"],"results/finalData/filtered_table.csv")
    output:
        os.path.join(config["general"]["output_dir"],"results/finalData/swarm_table.csv")
    conda:
        "../envs/merge_results.yaml"
    script:
        "../scripts/merge_swarm_results.py"
