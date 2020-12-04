import os

			
if config["classify"]["database"] == "PR2":
	rule download_pr2:
		output:
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "pr2db.fasta"),
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "pr2db.tax"),
			path = config["classify"]["db_path"]
		params:
			db_path=config["classify"]["db_path"],
			dnamol=config["classify"]["dnamol"],
			pr2_version=config["classify"]["pr2_version"]
		conda:
			"../envs/classify.yaml"
		shell:
			"""
			dir_name=$(dirname {params[0]});
			wget -N -P $dir_name/ https://github.com/pr2database/pr2database/releases/download/v{params[2]}/pr2_version_{params[2]}_{params[1]}_mothur.fasta.gz;
			wget -N -P $dir_name/ https://github.com/pr2database/pr2database/releases/download/v{params[2]}/pr2_version_{params[2]}_{params[1]}_mothur.tax.gz;
			gunzip -c $dir_name/pr2_version_{params[2]}_{params[1]}_mothur.fasta.gz > $dir_name/pr2db.fasta;	
			gunzip -c $dir_name/pr2_version_{params[2]}_{params[1]}_mothur.tax.gz > $dir_name/pr2db.tax;
			touch {output.path}
			"""

elif config["classify"]["database"] == "SILVA":
	rule download_silva_mothur:
		output:
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "silvadb.fasta"),
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "silvadb.tax"),
			path = config["classify"]["db_path"]
		params:
			db_path=config["classify"]["db_path"],
			silva_version=config["classify"]["silva_version"]
		conda:
			"../envs/classify.yaml"
		shell:
			"""
			dir_name=$(dirname {params[0]});
			wget -N -P $dir_name/ https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v{params[1]}.tgz;
			tar zxvf $dir_name/silva.nr_v{params[1]}.tgz -C $dir_name;
			touch {output.path};
			mothur "#filter.seqs(fasta=$dir_name/silva.nr_v{params[1]}.align, vertical = T)";
			mv $dir_name/silva.nr_v{params[1]}.filter.fasta $dir_name/silvadb.fasta;
			mv $dir_name/silva.nr_v{params[1]}.tax $dir_name/silvadb.tax 
			"""
			
elif config["classify"]["database"] == "RDP":
	rule download_rdp_mothur:
		output:
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "rdpdb.fasta"),
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "rdpdb.tax"),
			path = config["classify"]["db_path"]
		params:
			db_path=config["classify"]["db_path"]
			#rdp_version=config["classify"]["rdp_version"]
		conda:
			"../envs/classify.yaml"
		shell:
			"""
			dir_name=$(dirname {params[0]});
			wget -N -P $dir_name/ https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.rdp.tgz;
			tar zxvf $dir_name/trainset16_022016.rdp.tgz -C $dir_name;
			mv $dir_name/trainset16_022016.rdp/trainset16_022016.rdp.fasta $dir_name/rdpdb.fasta; mv $dir_name/trainset16_022016.rdp/trainset16_022016.rdp.tax $dir_name/rdpdb.tax;
			touch {output.path}
			"""

elif config["classify"]["database"] == "Unite":
	rule download_unite_mothur:
		output:
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "unitedb.fasta"),
			os.path.join(os.path.dirname(config["classify"]["db_path"]), "unitedb.tax"),
			path = config["classify"]["db_path"]
		params:
			db_path=config["classify"]["db_path"]
		conda:
			"../envs/classify.yaml"
		shell:
			"""
			dir_name=$(dirname {params[0]});
			wget -N -P $dir_name/ https://files.plutof.ut.ee/public/orig/CA/1C/CA1C2FA181FE7929BDBEE4644185280FFC332E31F18F721B241019BC43A31B64.gz;
			mv $dir_name/CA1C2FA181FE7929BDBEE4644185280FFC332E31F18F721B241019BC43A31B64.gz $dir_name/unite.gz;
			gunzip $dir_name/unite.gz; mv $dir_name/unite $dir_name/unite.tar;
			tar xvf $dir_name/unite.tar -C $dir_name; mv $dir_name/sh_mothur_release_s_all_04.02.2020/UNITEv8_sh_99_s_all.* $dir_name/;
			mv $dir_name/UNITEv8_sh_99_s_all.fasta $dir_name/unitedb.fasta; mv $dir_name/UNITEv8_sh_99_s_all.tax $dir_name/unitedb.tax;
			touch {output.path}
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
		mothur "#set.dir(blastdir=~/miniconda3/envs/natrix/bin);classify.seqs(fasta={input[0]}, template={input[1]}, taxonomy = {input[2]}, method={params.method}, ksize={params.ksize}, iters={params.iters}, search={params.search}, cutoff ={params.cutoff}, probs={params.probs}, numwanted={params.numwanted}, gapopen={params.gapopen}, gapextend={params.gapextend}, relabund={params.relabund}, output={params.output}, printlevel={params.printlevel})" -q > mothurout;
		mv results/finalData/representatives.*.{params[0]}.taxonomy results/finalData/classify_taxonomy.tsv;
		mv results/finalData/representatives.*.{params[0]}.tax.summary results/finalData/classify_summary.txt;
		mv mothur*logfile results/finalData/;
		rm mothurout;
		"""
