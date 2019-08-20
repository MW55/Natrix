rule all:
    input:
        "files/filtered_out_table.csv"

rule filtering:
    input:
        "files/unfiltered_table.hdf5"
    output:
        temp("files/filtered_table_temp.csv"),
        "files/filtered_out_table.csv"
    params:
        filter_method=config["merge"]["filter_method"],
        cutoff=config["merge"]["cutoff"]
    conda:
        "../envs/filtering.yaml"
    script:
        "../scripts/filtering.py"
