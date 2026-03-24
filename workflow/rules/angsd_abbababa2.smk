ABBABABA2_DIR = f"results/angsd_abbababa2/{output_prefix}"
ABBABABA2_SUMMARY_DIR = f"results/abbababa2/{output_prefix}"
ABBABABA2_FIG_DIR = f"figures/abbababa2/{output_prefix}"


rule prepare_abbababa2_inputs:
    """
    Assemble bamlist, sizeFile, and popNames for ABBABABA2 based on config populations and selected outgroups.
    """
    input:
        samples=config["samples"],
        sliced=(lambda wc: [] if not ABBABABA2_OUTGROUP_IDS else expand(f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam", sample_id=ABBABABA2_OUTGROUP_IDS))
    output:
        bamlist=f"results/bamlists/{output_prefix}/abbababa2/bamlist.txt",
        sizefile=f"{ABBABABA2_DIR}/sizeFile",
        popnames=f"{ABBABABA2_DIR}/popNames"
    params:
        bam_dir=config_bam_dir
    run:
        import os
        import pandas as pd

        if not ABBABABA2_OUTGROUP_IDS:
            raise ValueError("abbababa2.outgroup_samples is empty; at least one outgroup sample_id is required")

        df = pd.read_csv(input.samples, sep="\t")
        missing = [pop for pop in groups if df[df[group_col] == pop].empty]
        if missing:
            raise ValueError(f"No samples found for populations {missing} in {input.samples} for ABBABABA2")

        records = []
        for pop in groups:
            label = get_population_label(pop)
            subset = df[df[group_col] == pop]
            for sample in subset["sample"].astype(str):
                if sample in ABBABABA2_EXCLUDE_SAMPLES:
                    continue
                records.append((label, f"{params.bam_dir}/{sample}.bam"))

        for sid in ABBABABA2_OUTGROUP_IDS:
            label = resolve_abbababa2_outgroup_label(sid)
            records.append((label, f"{OUTGROUP_SLICED_DIR}/{sid}.bam"))
        os.makedirs(os.path.dirname(output.sizefile), exist_ok=True)

        pop_order = []
        counts = {}
        bam_paths = []
        for label, path in records:
            if label not in counts:
                counts[label] = 0
                pop_order.append(label)
            counts[label] += 1
            bam_paths.append(path)

        if len(pop_order) != len(set(pop_order)):
            raise ValueError("population_labels produce duplicate names for ABBABABA2. Ensure each population has a unique label.")

        write_bamlist(output.bamlist, bam_paths)

        with open(output.sizefile, "w") as handle:
            handle.write("\n".join(str(counts[label]) for label in pop_order) + "\n")

        with open(output.popnames, "w") as handle:
            handle.write("\n".join(pop_order) + "\n")


rule angsd_abbababa2:
    """
    Run ANGSD with -doAbbababa2 on intersecting sites.
    """
    input:
        bamlist=rules.prepare_abbababa2_inputs.output.bamlist,
        sizefile=rules.prepare_abbababa2_inputs.output.sizefile,
        popnames=rules.prepare_abbababa2_inputs.output.popnames,
        sites=f"results/intersect_sites/{output_prefix}/intersect.txt",
        scafs=f"results/intersect_sites/{output_prefix}/intersect.chr",
        sites_idx=f"results/intersect_sites/{output_prefix}/intersect.txt.bin"
    output:
        abbababa=f"{ABBABABA2_DIR}/bam.Angsd.abbababa2"
    log:
        f"logs/{output_prefix}/angsd_abbababa2.log"
    params:
        ref=REF,
        outprefix=f"{ABBABABA2_DIR}/bam.Angsd",
        extra=ABBABABA2_ANGSD_ARGS,
        minInd_ratio=get_minInd_ratio("abbababa2", get_minInd_ratio("global", None))
    threads: ABBABABA2_THREADS
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        MININD_OPT=""
        if [ -n "{params.minInd_ratio}" ] && [ "{params.minInd_ratio}" != "None" ]; then
          N=$(wc -l < {input.bamlist})
          r="{params.minInd_ratio}"
          MININD=$(awk -v n="$N" -v r="$r" 'BEGIN{{mi=int(n*r+0.5); if(mi<1) mi=1; print mi}}')
          MININD_OPT="-minInd $MININD"
        fi

        angsd -out {params.outprefix} -b {input.bamlist} \
              -ref {params.ref} -anc {params.ref} \
              -sites {input.sites} \
              -rf {input.scafs} \
              -sizeFile {input.sizefile} \
              {params.extra} $MININD_OPT \
              -nThreads {threads} \
              2> {log}
        """


rule abbababa2_estavgerror:
    input:
        abbababa=rules.angsd_abbababa2.output.abbababa,
        sizefile=rules.prepare_abbababa2_inputs.output.sizefile,
        popnames=rules.prepare_abbababa2_inputs.output.popnames
    output:
        observed=f"{ABBABABA2_DIR}/result.Observed.txt",
        transrem=f"{ABBABABA2_DIR}/result.TransRem.txt"
    params:
        angsd_prefix=f"{ABBABABA2_DIR}/bam.Angsd",
        out=f"{ABBABABA2_DIR}/result"
    log:
        f"logs/{output_prefix}/abbababa2_estAvgError.log"
    conda:
        "../envs/abbababa2_estavg.yaml"
    shell:
        """
        Rscript workflow/scripts/estAvgError.R \
            angsdFile="{params.angsd_prefix}" \
            out="{params.out}" \
            sizeFile="{input.sizefile}" \
            nameFile="{input.popnames}" \
            > {log} 2>&1
        """


rule summarize_abbababa2:
    input:
        observed=rules.abbababa2_estavgerror.output.observed,
        transrem=rules.abbababa2_estavgerror.output.transrem
    output:
        summary=f"{ABBABABA2_SUMMARY_DIR}/df.dstat.csv"
    conda:
        "../envs/python_pandas.yaml"
    script:
        "../scripts/summarize_abbababa2.py"


rule plot_abbababa2:
    input:
        summary=rules.summarize_abbababa2.output.summary
    output:
        pdf=f"{ABBABABA2_FIG_DIR}/abbababa2.pdf"
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot_abbababa2.R"
