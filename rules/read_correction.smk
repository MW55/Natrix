if config['dataset']['nanopore']:
    ## fastq to fasta
    rule fastq2fasta:
        input:
            expand(os.path.join(config["general"]["output_dir"],"pychopper_merged/{{sample}}_{{unit}}_R{read}.fastq"),read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
        conda:
            "../envs/read_correction.yaml"
        shell:
            "fastq_to_fasta -i {input} -o {output}"
    ## cdhit
    rule cd_hit:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads)
        conda: "../envs/read_correction.yaml"
        shell:
            " cd-hit-est -i {input} -o {output} -c 0.9 -d 0 -M 0 -T 20"

    ## alighnment with minimap on reads
    rule minimap_align:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads)
        conda:
             "../envs/read_correction.yaml"
        shell:
            "minimap2 -ax map-ont -t 20 {input[0]}  {input[1]} > {output}"

    ## polishing with racon
    rule racon_polishing:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align.sam"), read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.sam"),read=reads)
        conda:
            "../envs/read_correction.yaml"
        shell:
            "racon  -u -f -w 50 -q 9  -t 20 {input[0]} {input[1]} {input[2]}  > polished.fasta"

    ## medaka polishing
    rule medaka_polishing:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.sam"), read=reads)
        params:
            prefix="consensus",
            out_dir=expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads)

        conda: "../envs/medaka.yaml"
        shell:
            """
            medaka_consensus -i {input[0]} -d {input[1]} -o {params.out_dir} -t 20"
            """
