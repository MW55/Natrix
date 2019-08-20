rule unfiltered_table:
    input:
        expand("results/finalData/{unit.sample}_{unit.unit}.nonchimera.fasta", unit=units.reset_index().itertuples())
    output:
        "results/finalData/unfiltered_table.csv",
        temp("results/finalData/unfiltered_table.hdf5")
    conda:
        "../envs/unfiltered_table.yaml"
    script:
        "../scripts/unfiltered_table.py"

rule filtering:
    input:
        "results/finalData/unfiltered_table.hdf5"
    output:
        temp("results/finalData/filtered_table_temp.csv"),
        "results/finalData/filtered_out_table.csv"
    params:
        filter_method=config["merge"]["filter_method"],
        cutoff=config["merge"]["cutoff"]
    conda:
        "../envs/filtering.yaml"
    script:
        "../scripts/filtering.py"

rule ampliconduo:
    input:
        filtered_table="results/finalData/filtered_table_temp.csv",
        unfiltered_table="results/finalData/unfiltered_table.csv"
    output:
        "results/finalData/figures/AmpliconDuo.RData"
    params:
        plot_ampduo=config["merge"]["plot_AmpDuo"],
        saving_format=config["merge"]["save_format"],
        p_corr=config["merge"]["ampli_corr"]
    conda:
        "../envs/ampliconduo.yaml"
    script:
        "../scripts/ampliconduo.R"
