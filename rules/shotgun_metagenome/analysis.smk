rule diversity_analysis:
    input:
        taxonomy_matrix="results/taxonomy/taxonomy_abundance_matrix.tsv"
    output:
        alpha_diversity="results/diversity/alpha_diversity.tsv",
        beta_diversity="results/diversity/beta_diversity_pcoa.tsv"
    conda:
        "../../envs/diversity_analysis.yaml"
    script:
        "../../scripts/shotgun/compute_diversity.py"

rule plot_taxonomy_functional:
    input:
        taxonomy_matrix="results/taxonomy/taxonomy_abundance_matrix.tsv",
        functional_matrix="results/functional/functional_abundance_matrix.tsv"
    output:
        taxonomy_plot="results/visualization/taxonomy_barplot.png",
        functional_heatmap="results/visualization/functional_heatmap.png"
    conda:
        "../../envs/visualization.yaml"
    script:
        "../../scripts/shotgun/plot_taxonomy_functional.py"

rule plot_beta_diversity:
    input:
        beta_diversity="results/diversity/beta_diversity_pcoa.tsv"
    output:
        pcoa_plot="results/visualization/beta_diversity_pcoa.png"
    conda:
        "../../envs/visualization.yaml"
    script:
        "../../scripts/shotgun/plot_beta_diversity.py"