# ngsadmix.smk

rule run_ngsadmix:
    input:
        beagle = rules.angsd_global_unrelated_unlinked.output.beagle
    output:
        expand(f"results/admixture/{output_prefix}/K{{k}}_{{r}}.qopt", 
               k=range(1, config["ngsadmix"]["maxK"] + 1), 
               r=range(1, config["ngsadmix"]["n_replicates"] + 1)),
               expand(f"results/admixture/{output_prefix}/K{{k}}_{{r}}.log", 
               k=range(1, config["ngsadmix"]["maxK"] + 1), 
               r=range(1, config["ngsadmix"]["n_replicates"] + 1))
    params:
        outdir = f"results/admixture/{output_prefix}",
        threads = config["ngsadmix"]["threads"],
        maxiter = config["ngsadmix"]["maxiter"],
        minmaf = config["ngsadmix"]["minmaf"],
        maxK = config["ngsadmix"]["maxK"],
        reps = config["ngsadmix"]["n_replicates"]
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        for K in $(seq 1 {params.maxK}); do
          for REP in $(seq 1 {params.reps}); do
            NGSadmix -likes {input.beagle} \
                     -K $K \
                     -outfiles {params.outdir}/K${{K}}_${{REP}} \
                     -P {params.threads} \
                     -maxiter {params.maxiter} \
                     -minMaf {params.minmaf} \
                     -printInfo 1
          done
        done
        """

rule samples_with_outgroups_admix:
    input:
        samples=config["samples"],
        outgroups=config["outgroups"]
    output:
        samples_aug=f"results/metadata/{output_prefix}/samples_plus_outgroups.tsv"
    run:
        import pandas as pd, os
        os.makedirs(os.path.dirname(output.samples_aug), exist_ok=True)
        s = pd.read_csv(input.samples, sep="\t")
        og = pd.read_csv(input.outgroups, sep="\t")
        # build rows with columns matching samples: assume 'sample' exists; fill group_col with 'taxon'
        cols = s.columns.tolist()
        if "sample" not in cols:
            raise ValueError("samples TSV must include a 'sample' column")
        # prepare new rows
        rows = []
        for _, r in og.iterrows():
            row = {c: None for c in cols}
            row["sample"] = str(r["sample_id"]) if "sample_id" in og.columns else None
            if group_col in cols and "taxon" in og.columns:
                row[group_col] = str(r["taxon"]) if pd.notna(r["taxon"]) else None
            rows.append(row)
        aug = pd.concat([s, pd.DataFrame(rows)], ignore_index=True)
        aug.to_csv(output.samples_aug, sep="\t", index=False)

rule plot_admixture:
    input:
        qopt=expand(
            f"results/admixture/{output_prefix}/K{{k}}_{{r}}.qopt",
            k=range(1, config["ngsadmix"]["maxK"]+1),
            r=range(1, config["ngsadmix"]["n_replicates"]+1)
        ),
        logs=expand(
            f"results/admixture/{output_prefix}/K{{k}}_{{r}}.log",
            k=range(1, config["ngsadmix"]["maxK"]+1),
            r=range(1, config["ngsadmix"]["n_replicates"]+1)
        ),
        bamlist=rules.make_bamlist_unrelated_analysis.output.bamlist,
        samples=rules.samples_with_outgroups_admix.output.samples_aug,
        beagle=rules.angsd_global_unrelated_unlinked.output.beagle
    output:
        delta="results/admixture/{output_prefix}/deltaK.pdf",
        plots="results/admixture/{output_prefix}/admixture_plots.pdf"
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot_admixture.R"
