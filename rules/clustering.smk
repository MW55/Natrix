# dada runs per sample after preprocessing and before chimera removal and AmpliconDuo
rule DADA2:
    input:
        fwd = expand(
            "results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_1_cut.fastq", unit=units.reset_index().itertuples()),
        rev = expand(
            "results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_2_cut.fastq", unit=units.reset_index().itertuples()) if config["merge"]["paired_End"] == True else []
    output:
        expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_dada.fasta", unit=units.reset_index().itertuples())
    params:
        paired_end=config["merge"]["paired_End"],
        minoverlap=config["qc"]["minoverlap"],
        splitsamples=config["merge"]["filter_method"]
    conda:
        "../envs/dada2.yaml"
    log:
        "results/logs/dada2.log"
    script:
        "../scripts/dada2.R"

# SWARM clustering runs on all samples after AmpliconDuo
rule write_fasta:
    input:
        "results/finalData/filtered_table_temp.csv"
    output:
        "results/finalData/filtered.fasta",
        "results/finalData/filtered_table.csv"
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

if config["general"]["sequencing"] == "Nanopore":
    rule vsearch_otu:
        input:
            "results/finalData/filtered.fasta"
        output:
            "results/finalData/representatives.fasta",
            "results/finalData/merged.uc"
        params:
            cutoff=config["merge"]["vsearch_clust_id"],
            output=config["merge"]["vearch_clust_output"]
        threads: config["general"]["cores"]
        conda:
            "../envs/vsearch.yaml"
        shell:
            "vsearch --cluster_fast {input} --{params.output} {output[0]} --id {params.cutoff} --uc {output[1]} --sizein --clusterout_id --threads {threads}"
    
    rule vsearch_otu_results:
        input:
            merged="results/finalData/merged.uc",
            consensus="results/finalData/representatives.fasta",
            final_table_path2="results/finalData/filtered_table.csv"
        output:
            swarm_table="results/finalData/swarm_table.csv",
            all_out="results/finalData/swarm_table_all.csv",
            consensus_filtered="results/finalData/representatives_filtered.fasta"
        params:
            min_sequences=config["merge"]["vsearch_clust_min"],
            consensus=config["merge"]["vsearch_clust_output"]
        conda:
            "../envs/merge_results.yaml"
        script:
            "../scripts/merge_cluster_results.py"
        

else:
    rule swarm:
        input:
            "results/finalData/filtered.fasta"
        output:
            "results/finalData/representatives.fasta",
            temp("results/finalData/merged.swarms")
        threads: config["general"]["cores"]
        conda:
            "../envs/swarm.yaml"
        shell:
            "swarm -t {threads} -d 1 -f -z -w {output[0]} < {input} > {output[1]}"

    rule swarm_results:
        input:
            merged="results/finalData/merged.swarms",
            final_table_path2="results/finalData/filtered_table.csv"
        output:
            "results/finalData/swarm_table.csv"
        conda:
            "../envs/merge_results.yaml"
        script:
            "../scripts/merge_swarm_results.py"
