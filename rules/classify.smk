import os

rule download_pr2:
	output:
		os.path.join(os.path.dirname(config["classify"]["db_path"]), "pr2.db.fasta"),
		os.path.join(os.path.dirname(config["classify"]["db_path"]), "pr2.db.tax")
	params:
		db_path=config["classify"]["db_path"],
		dnamol=config["classify"]["dnamol"]
	conda:
		"../envs/classify.yaml"
	shell:
		"""
		dir_name=$(dirname {params[0]});
		wget -N -P $dir_name/ https://github.com/pr2database/pr2database/releases/download/v4.12.0/pr2_version_4.12.0_{params[1]}_mothur.fasta.gz;
		wget -N -P $dir_name/ https://github.com/pr2database/pr2database/releases/download/v4.12.0/pr2_version_4.12.0_{params[1]}_mothur.tax.gz;
		gunzip -c $dir_name/pr2_version_4.12.0_{params[1]}_mothur.fasta.gz > $dir_name/pr2.db.fasta;	
		gunzip -c $dir_name/pr2_version_4.12.0_{params[1]}_mothur.tax.gz > $dir_name/pr2.db.tax;
		"""
		
rule classify:
	input:
		"results/finalData/representatives.fasta" if config["general"]["seq_rep"] == "OTU" else "results/finalData/filtered.fasta", 
		expand(config["classify"]["db_path"] + "{file_extension}", file_extension=[".fasta", ".tax"])
	output:
		"results/finalData/classify_taxonomy.tsv",
		"results/finalData/classify_summary.txt"
	params: 
		method=config["classify"]["method"],
		ksize=config["classify"]["ksize"],
		iters=config["classify"]["iters"],
		cutoff=config["classify"]["cutoff"],
		probs=config["classify"]["probs"],
		search=config["classify"]["search"],
		numwanted=config["classify"]["numwanted"],
		gapopen=config["classify"]["gapopen"],
		gapextend=config["classify"]["gapextend"],
		relabund=config["classify"]["relabund"],
		output=config["classify"]["output"],
		printlevel=config["classify"]["printlevel"]
	threads: 
		config["general"]["cores"]
	conda:
		"../envs/classify.yaml"
	
	shell:
		"""
		mothur "#classify.seqs(fasta={input[0]}, template={input[1]}, taxonomy = {input[2]}, method={params[0]}, ksize={params[1]}, cutoff ={params[2]}, gapopen={params[7]})" > mothurlog; #mothur creates another log file
		mv results/finalData/representatives.db.wang.taxonomy results/finalData/classify_taxonomy.tsv;
		mv results/finalData/representatives.db.wang.tax.summary results/finalData/classify_summary.txt;
		mv mothur*logfile results/finalData/;
		rm mothurlog;
		"""
