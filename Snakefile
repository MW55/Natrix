configfile: "test_data.yaml"

fastq_dir = config["general"]["filename"]
IDS, = glob_wildcards(fastq_dir + "/{id}.fastq.gz")

rule fastqc:
	input:
		fastq = expand(fastq_dir + "/{id}.fastq.gz", id=IDS)
	output:
		"results/qc/ {input.id} _fastqc.html",
		"results/qc/ {input.id} _fastqc.zip" 
	shell:
		"fastqc -o results/qc/ {input.fastq}"
	
	#output:
		#"results/qc/ {input.fastq}"

