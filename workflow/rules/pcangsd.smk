# rules/pcangsd.smk
rule run_pcangsd:
    input:
        beagle = rules.angsd_global_unrelated_unlinked.output.beagle
    output:
        pca = f"results/pca/{output_prefix}/pcangsd.cov",
        pca_tree = f"results/pca/{output_prefix}/pcangsd.tree.cov"
    log:
        f"logs/{output_prefix}/pcangsd.log"
    params:
        outprefix = f"results/pca/{output_prefix}/pcangsd",
        iter = config["pcangsd"]["iter"],
        minmaf = config["pcangsd"]["minmaf"]
    singularity:
        f"{config_singularity_dir}/pcangsd_1.35.sif"
    shell:
        """
        pcangsd --beagle {input.beagle} \
                -o {params.outprefix} \
                --iter {params.iter} \
                --maf {params.minmaf} \
                --admix --tree --selection --snp_weights --sites_save 2> {log}
        """

rule pca_selection:
    input:
        selection=f"results/pca/{output_prefix}/pcangsd.selection",
        pca_sites=f"results/pca/{output_prefix}/pcangsd.sites"
    output:
        selection_csv = f"results/pca/{output_prefix}/selection.csv",
        qqplot = f"results/pca/{output_prefix}/pca_{output_prefix}_qq.png"
    log:
        f"logs/{output_prefix}/pca_selection.log"
    conda:
        "../envs/pca_selection.yamlß"
    shell:
        """
        Rscript workflow/scripts/run_selection.R \
            {input.selection} {input.pca_sites}\
            results/pca/{output_prefix}/pcangsd.sites \
            {output.selection_csv} {output.qqplot}
        """

rule plot_pca:
    input:
        cov=rules.run_pcangsd.output.pca,
        eig=rules.run_pcangsd.output.pca_tree,
        bamlist=rules.make_bamlist_unrelated_analysis.output.bamlist,
        samples=f"results/metadata/{output_prefix}/samples_plus_outgroups.tsv"
    output:
        pdf=f"results/pca/{output_prefix}/pca_plots.pdf"
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot_pca.R"


## samples_with_outgroups_pca removed; use the producer in ngsadmix.smk
