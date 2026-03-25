# rules/pcangsd.smk

PCA_DIR = f"results/pca/{output_prefix}"
PCA_FIG_DIR = f"figures/exploratory/pca/{output_prefix}"
PCA_SELECTION_FIG_DIR = f"{PCA_FIG_DIR}/selection"

rule run_pcangsd:
    """
    Run PCAngsd on the unrelated, LD-pruned genotype-likelihood matrix.
    """
    input:
        beagle = rules.angsd_global_unrelated_unlinked.output.beagle
    output:
        pca = f"{PCA_DIR}/pcangsd.cov",
        pca_tree = f"{PCA_DIR}/pcangsd.tree.cov",
        selection = f"{PCA_DIR}/pcangsd.selection",
        sites = f"{PCA_DIR}/pcangsd.sites"
    log:
        f"logs/{output_prefix}/pcangsd.log"
    params:
        outprefix = f"{PCA_DIR}/pcangsd",
        iter = config["pcangsd"]["iter"],
        minmaf = config["pcangsd"]["minmaf"]
    threads: PCANGSD_THREADS
    singularity:
        f"{config_singularity_dir}/pcangsd_1.35.sif"
    shell:
        """
        pcangsd --beagle {input.beagle} \
                -o {params.outprefix} \
                --iter {params.iter} \
                --threads {threads} \
                --maf {params.minmaf} \
                --admix --tree --selection --snp_weights --sites_save 2> {log}
        """

rule pca_selection:
    """
    Summarize PCAngsd site-selection output and render the QQ diagnostic plot.
    """
    input:
        selection=rules.run_pcangsd.output.selection,
        pca_sites=rules.run_pcangsd.output.sites
    output:
        selection_csv = f"{PCA_DIR}/selection.csv",
        qqplot = f"{PCA_SELECTION_FIG_DIR}/pca_{output_prefix}_qq.png"
    log:
        f"logs/{output_prefix}/pca_selection.log"
    conda:
        "../envs/r_plotting.yaml"
    shell:
        """
        Rscript workflow/scripts/run_selection.R \
            {input.selection} {input.pca_sites}\
            {output.selection_csv} {output.qqplot}
        """

rule plot_pca:
    """
    Plot PCA coordinates for the configured axes and write a multi-page PDF.
    """
    input:
        cov=rules.run_pcangsd.output.pca,
        eig=rules.run_pcangsd.output.pca_tree,
        bamlist=rules.make_bamlist_unrelated_analysis.output.bamlist,
        samples=f"results/metadata/{output_prefix}/samples_plus_outgroups.tsv"
    output:
        pdf=f"{PCA_FIG_DIR}/pca_plots.pdf"
    conda:
        "../envs/r_plotting.yaml"
    script:
        "../scripts/plot_pca.R"


## samples_with_outgroups_pca removed; use the producer in ngsadmix.smk

PCANGSD_TARGETS = [
    rules.run_pcangsd.output.pca,
    rules.run_pcangsd.output.pca_tree,
    rules.run_pcangsd.output.selection,
    rules.run_pcangsd.output.sites,
    rules.pca_selection.output.selection_csv,
    rules.pca_selection.output.qqplot,
    rules.plot_pca.output.pdf,
] if PCANGSD_ENABLED else []
