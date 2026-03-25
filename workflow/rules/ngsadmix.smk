"""
NGSadmix clustering and plotting workflow on unrelated, LD-pruned data.

Inputs:
- Unrelated, unlinked BEAGLE from angsd_global.smk
- Sample and outgroup metadata for plot annotation

Config keys used:
- ngsadmix.enabled
- ngsadmix.maxK
- ngsadmix.n_replicates
- ngsadmix.maxiter
- ngsadmix.minmaf
- populations
- group_col
- threads.analyses.ngsadmix
"""

ADMIX_DIR = f"results/admixture/{output_prefix}"
ADMIX_FIG_DIR = f"figures/exploratory/admixture/{output_prefix}"
ADMIX_MAX_K = int(config["ngsadmix"]["maxK"])
ADMIX_REPLICATES = int(config["ngsadmix"]["n_replicates"])

rule run_ngsadmix:
    """
    Run NGSadmix across the configured K range and replicate count.
    """
    input:
        beagle = rules.angsd_global_unrelated_unlinked.output.beagle
    output:
        expand(f"{ADMIX_DIR}/K{{k}}_{{r}}.qopt", k=range(1, ADMIX_MAX_K + 1), r=range(1, ADMIX_REPLICATES + 1)),
        expand(f"{ADMIX_DIR}/K{{k}}_{{r}}.log", k=range(1, ADMIX_MAX_K + 1), r=range(1, ADMIX_REPLICATES + 1))
    params:
        outdir = ADMIX_DIR,
        maxiter = config["ngsadmix"]["maxiter"],
        minmaf = config["ngsadmix"]["minmaf"],
        maxK = ADMIX_MAX_K,
        reps = ADMIX_REPLICATES
    threads: NGSADMIX_THREADS
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        for K in $(seq 1 {params.maxK}); do
          for REP in $(seq 1 {params.reps}); do
            NGSadmix -likes {input.beagle} \
                     -K $K \
                     -outfiles {params.outdir}/K${{K}}_${{REP}} \
                     -P {threads} \
                     -maxiter {params.maxiter} \
                     -minMaf {params.minmaf} \
                     -printInfo 1
          done
        done
        """

rule plot_admixture:
    """
    Plot DeltaK and admixture barplots for the best replicate at each K.
    """
    input:
        qopt=expand(
            f"{ADMIX_DIR}/K{{k}}_{{r}}.qopt",
            k=range(1, ADMIX_MAX_K + 1),
            r=range(1, ADMIX_REPLICATES + 1)
        ),
        logs=expand(
            f"{ADMIX_DIR}/K{{k}}_{{r}}.log",
            k=range(1, ADMIX_MAX_K + 1),
            r=range(1, ADMIX_REPLICATES + 1)
        ),
        bamlist=rules.make_bamlist_unrelated_analysis.output.bamlist,
        samples=rules.samples_with_outgroups_metadata.output.samples_aug,
        beagle=rules.angsd_global_unrelated_unlinked.output.beagle
    output:
        delta=f"{ADMIX_FIG_DIR}/deltaK.pdf",
        plots=f"{ADMIX_FIG_DIR}/admixture_plots.pdf"
    conda:
        "../envs/r_plotting.yaml"
    script:
        "../scripts/plot_admixture.R"


NGSADMIX_TARGETS = [
    *list(rules.run_ngsadmix.output),
    rules.samples_with_outgroups_metadata.output.samples_aug,
    rules.plot_admixture.output.delta,
    rules.plot_admixture.output.plots,
] if NGSADMIX_ENABLED else []
