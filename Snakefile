import pandas as pd

configfile: "test_data.yaml"
fastq_dir = config["general"]["filename"]
samples = pd.read_table(config["general"]["samples"], index_col="sample")
units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],	
	dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])


rule fastqc:
	input:
		sample = expand("test_data/{sample}_{unit}_R{group}.fastq.gz", group = [1, 2])

	output:
		"results/qc/{sample}_{unit}_R{group}_fastqc.html",
		"results/qc/{sample}_{unit}_R{group}_fastqc.zip" 
	shell:
		"fastqc -o results/qc/ {input.fastq}"
