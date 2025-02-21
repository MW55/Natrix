rule aggregate_taxonomy_results:
    """
    Aggregates multiple Bracken reports into one table for downstream diversity analysis.
    Each input file corresponds to one sample.
    """
    input:
        [f"results/taxonomy/{row.sample}_{row.unit}_bracken_report.txt"
            for row in units.reset_index().itertuples()]
    output:
        taxonomy_matrix="results/taxonomy/taxonomy_abundance_matrix.tsv"
    conda:
        "../../envs/aggregation.yaml"
    script:
        "../../scripts/shotgun/aggregate_bracken_reports.py"


rule aggregate_functional_results:
    input:
        gene_families=expand("results/functional/{sample}_{unit}_genefamilies.tsv", sample=list(units.reset_index()["sample"]), unit=list(units.reset_index()["unit"])),
        pathways=expand("results/functional/{sample}_{unit}_pathways.tsv",sample=list(units.reset_index()["sample"]), unit=list(units.reset_index()["unit"]))
    output:
        functional_matrix="results/functional/functional_abundance_matrix.tsv"
    conda:
        "../../envs/aggregation.yaml"
    shell:
        """
        humann_join_tables --input results/functional/ --output {output.functional_matrix} --file_name genefamilies.tsv
        humann_join_tables --input results/functional/ --output {output.functional_matrix} --file_name pathways.tsv
        """
