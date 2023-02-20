import os
if config['dataset']['nanopore']:
    rule pychop:
        input:
            fastq=expand(os.path.join(config["general"]["filename"],"{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
        output:
            out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/output/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            pdf_out=expand(os.path.join(config["general"]["output_dir"],"pychopper/reports/{{sample}}_{{unit}}_R{read}.pdf"), read=reads),
            unclass_fastq= expand(os.path.join(config["general"]["output_dir"],"pychopper/unclassified/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            rescue_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/rescued/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
        params:
            primer = config["general"]["primertable"],
            qual =config["pychop"]["qual"]
        conda:
            "../envs/pychopper.yaml"
        shell:
                """
                forward_primer=$(sed -n "2p" {params.primer} | cut -d, -f4);
                reverse_primer=$(sed -n "2p" {params.primer} | cut -d, -f7);
                echo "">SSP\\n$forward_primer\\n>VNP\\n$reverse_primer"" > custom_primers.fasta;
                echo "+:SSP,-VNP|-:VNP,-SSP" > config_primers.txt ;
                pychopper -m edlib -b custom_primers.fasta -Q {params.qual}  -c config_primers.txt -r {output.pdf_out} -u {output.unclass_fastq} -w {output.rescue_fastq} {input.fastq} {output.out_fastq};
                """

    rule pychopper_rescue:
        input:
            unclass_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper/unclassified/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
        output:
            unclass_out_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper_unclass/output/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            unclass_unclass_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper_unclass/unclassified/{{sample}}_{{unit}}_R{read}.fastq"),  read=reads),
            unclass_rescue_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper_unclass/rescued/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            unclass_pdf=expand(os.path.join(config["general"]["output_dir"],"pychopper_unclass/reports/{{sample}}_{{unit}}_R{read}.pdf"),  read=reads)
        params:
            primer=config["general"]["primertable"],
            qual =config["pychop"]["qual"]
        conda:
            "../envs/pychopper.yaml"
        shell:
            """
            forward_primer=$(sed -n "2p" {params.primer} | cut -d, -f4);
            reverse_primer=$(sed -n "2p" {params.primer} | cut -d, -f7);
            echo "">SSP\\n$forward_primer\\n>VNP\\n$reverse_primer"" > custom_primers.fasta;
            echo "+:SSP,-VNP|-:VNP,-SSP" > config_primers.txt;
    
            pychopper \
            -m edlib \
            -x rescue \
            -b custom_primers.fasta \
            -Q {params.qual}
            -c config_primers.txt \
            -r {output.unclass_pdf} \
            -u {output.unclass_unclass_fastq} \
            -w {output.unclass_rescue_fastq} \
            {input.unclass_fastq} \
            {output.unclass_out_fastq}
            """
    rule merge_pychopper:
        input:
            out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/output/{{sample}}_{{unit}}_R{read}.fastq"),  read=reads),
            unclass_out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper_unclass/output/{{sample}}_{{unit}}_R{read}.fastq"),  read=reads)
        output:
            pychopper_merged = expand(os.path.join(config["general"]["output_dir"], "pychopper_merged/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
        shell:
            """
            cat {input.out_fastq} {input.unclass_out_fastq} > {output.pychopper_merged}
            """



