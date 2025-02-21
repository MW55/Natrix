rule download_host_index:
    output:
        expand("database/host/human_index.{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    params:
        index_prefix="database/host/human_index",
        fasta_url="ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p database/host
        if [ ! -f {params.index_prefix}.1.bt2 ]; then
            wget -O database/host/human_reference.fa.gz {params.fasta_url}
            gunzip database/host/human_reference.fa.gz
            bowtie2-build database/host/human_reference.fa {params.index_prefix}
        fi
        """



rule host_removal:
    input:
        fq1="results/assembly/{sample}_{unit}/{sample}_{unit}_1_cut.fastq",
        fq2="results/assembly/{sample}_{unit}/{sample}_{unit}_2_cut.fastq",
        index="database/host/human_index.1.bt2"
    output:
        fq1="results/filtered/{sample}_{unit}_clean.1",
        fq2="results/filtered/{sample}_{unit}_clean.2"
    params:
        bowtie2_params=config["shotgun"]["host_removal"]["bowtie2_params"],
        threads=config["shotgun"]["host_removal"]["threads"],
        fq_prefix="results/filtered/{sample}_{unit}_clean"
    conda:
        "../../envs/host_removal.yaml"
    shell:
        """
        bowtie2 --quiet -x database/host/human_index -1 {input.fq1} -2 {input.fq2} \
        --un-conc {params.fq_prefix} \
        {params.bowtie2_params} --threads {params.threads}
        """