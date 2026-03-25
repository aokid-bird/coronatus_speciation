# rules/angsd_raxml.smk

RAXML_DIR = f"results/raxml/{output_prefix}"
RAXML_FIG_DIR = f"figures/exploratory/raxml/{output_prefix}"
RAXML_OUTGROUP_SPECIES = RAXML_CFG.get("outgroup_species")

rule make_bamlist_raxml_downsampled:
    """
    Downsample unrelated ingroup BAMs per population before RAxML analyses.
    """
    input:
        ingroup_unrel = rules.make_bamlist_unrelated.output.bamlist,
        samples=config["samples"]
    output:
        bamlist=f"results/bamlists/{output_prefix}/raxml_downsampled/bamlist.txt"
    params:
        group_col=group_col,
        populations=groups,
        max_per_pop=ANGSD_RAXML_DOWNSAMPLE_MAX,
        exclude=ANGSD_RAXML_DOWNSAMPLE_EXCLUDE,
        use_all=ANGSD_RAXML_DOWNSAMPLE_USE_ALL,
        seed=ANGSD_RAXML_DOWNSAMPLE_SEED
    run:
        import os
        import pandas as pd
        from pathlib import Path
        from random import Random

        meta = pd.read_csv(input.samples, sep="\t")
        if "sample" not in meta.columns:
            raise ValueError("'sample' column not found in samples metadata.")
        if params.group_col not in meta.columns:
            raise ValueError(
                f"Group column '{params.group_col}' missing from samples metadata."
            )

        sample_to_group = (
            meta.set_index("sample")[params.group_col]
            .astype(str)
            .to_dict()
        )

        exclude = {str(s) for s in params.exclude}

        bam_paths = pd.read_csv(input.ingroup_unrel, header=None)[0].tolist()
        entries = []  # (order, sample_id, group, path)
        for order, path in enumerate(bam_paths):
            sample_id = os.path.splitext(os.path.basename(str(path)))[0]
            if sample_id in exclude:
                continue
            group = sample_to_group.get(sample_id)
            if group is None:
                raise ValueError(
                    f"Sample '{sample_id}' from bamlist missing in samples metadata."
                )
            if params.populations and group not in params.populations:
                continue
            entries.append((order, sample_id, group, str(path)))

        grouped = {}
        for entry in entries:
            grouped.setdefault(entry[2], []).append(entry)

        def resolve_max(pop):
            max_cfg = params.max_per_pop
            if params.use_all:
                return None
            if max_cfg is None:
                return None
            if isinstance(max_cfg, dict):
                if pop in max_cfg and max_cfg[pop] is not None:
                    return int(max_cfg[pop])
                if "default" in max_cfg and max_cfg["default"] is not None:
                    return int(max_cfg["default"])
                return None
            return int(max_cfg)

        rng = Random(params.seed if params.seed is not None else 0)

        selected = []
        for pop in params.populations:
            candidates = grouped.get(pop, [])
            if not candidates:
                continue
            limit = resolve_max(pop)
            if limit is None or limit >= len(candidates):
                chosen = candidates
            elif limit <= 0:
                chosen = []
            else:
                idxs = rng.sample(range(len(candidates)), k=limit)
                chosen = [candidates[i] for i in sorted(idxs)]
            selected.extend(chosen)

        if params.use_all and not selected:
            selected = entries

        selected_sorted = sorted(selected, key=lambda x: x[0])
        if not selected_sorted:
            raise ValueError(
                "No samples selected for RAxML; check downsampling configuration or exclusions."
            )

        out_dir = Path(output.bamlist).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.Series([entry[3] for entry in selected_sorted]).to_csv(
            output.bamlist, index=False, header=False
        )

rule make_bamlist_raxml_analysis:
    """
    Bamlist for RAxML ANGSD; start from downsampled ingroup and optionally append sliced outgroups.
    """
    input:
        ingroup_downsampled = rules.make_bamlist_raxml_downsampled.output.bamlist,
        # Ensure sliced outgroup BAMs exist when including outgroups
        sliced=(lambda wc: sliced_outgroup_inputs(ANGSD_RAXML_OUTGROUP_IDS))
    output:
        bamlist = f"results/bamlists/{output_prefix}/raxml_analysis/bamlist.txt"
    run:
        import pandas as pd
        ing = pd.read_csv(input.ingroup_downsampled, header=None)[0].tolist()
        write_bamlist(output.bamlist, list(ing) + outgroup_bam_paths(ANGSD_RAXML_OUTGROUP_IDS))


rule angsd_raxml:
    """
    Run ANGSD on the RAxML analysis sample set to produce genotype likelihoods and BCF.
    """
    input:
        bamlist=rules.make_bamlist_raxml_analysis.output.bamlist,
        sites=f"results/intersect_sites/{output_prefix}/intersect.txt",
        scafs=f"results/intersect_sites/{output_prefix}/intersect.chr",
        sites_idx = f"results/intersect_sites/{output_prefix}/intersect.txt.bin"
    output:
        geno=f"results/angsd_raxml/{output_prefix}/gl.geno.gz",
        bcf=f"results/angsd_raxml/{output_prefix}/gl.bcf"
    log:
        f"logs/{output_prefix}/angsd_raxml.log"
    params:
        ref=REF,
        outprefix=f"results/angsd_raxml/{output_prefix}/gl",
        extra=config["angsd_common_args"].strip() + " " + config["angsd_args"]["raxml"].strip(),
        minInd_ratio=get_minInd_ratio("raxml", get_minInd_ratio("global", None))
    threads: ANGSD_RAXML_THREADS
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        # Compute optional -minInd based on ratio and bamlist size
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
              {params.extra} $MININD_OPT \
              -nThreads {threads} \
              2> {log}
        """

rule catg_format:
    """
    Convert the ANGSD VCF output into CATG format for RAxML-ng.
    """
    input:
        vcf=f"results/angsd_raxml/{output_prefix}/gl.vcf.gz",
        bamlist=rules.make_bamlist_raxml_analysis.output.bamlist
    output:
        raxmlcatg=f"{RAXML_DIR}/rxmlcatg.txt"
    params:
        maxSize = 4 * 1024**3 # max size for future.apply
    threads:config['raxml']['threads_catg']
    resources:
        mem_mb=200000,
        runtime=2880
    conda:
        "../envs/vcfR.yaml"
    script:
        "../scripts/catg_formatter.R"

rule raxml_ng:
    """
    Run RAxML-ng with bootstrap support on the CATG alignment.
    """
    input:
        raxmlcatg=rules.catg_format.output.raxmlcatg,
        # also pass bamlist to derive exact CATG labels for outgroup option
        bamlist=rules.make_bamlist_raxml_analysis.output.bamlist
    output:
        raxout=f"{RAXML_DIR}/rxmlcatg.txt.raxml.bootstraps",
        raxsup=f"{RAXML_DIR}/rxmlcatg.txt.raxml.support"
    threads:config['raxml']['threads_run']
    resources:
        mem_mb=200000,
        runtime=4320
    params:
        model=config['raxml']['model'],
        bs=config['raxml']['bs'],
        # Pre-render optional outgroup argument from actual bamlist labels
        outgroup_opt=(
            (lambda wildcards, input: outgroup_option_from_bamlist(
                input.bamlist,
                species=RAXML_OUTGROUP_SPECIES,
                sample_ids=ANGSD_RAXML_OUTGROUP_IDS,
             )) if ANGSD_RAXML_OUTGROUP_IDS else ""
        )
    conda:
        "../envs/raxml_ng.yaml"
    shell:
        """
        raxml-ng --msa {input.raxmlcatg} \
            --msa-format CATG --prob-msa on \
            --all --model {params.model} \
            --bs-trees {params.bs} \
            {params.outgroup_opt} \
            --threads {threads}
        """

rule plot_raxml:
    """
    Plot RAxML-ng bootstrap tree
    """
    input:
        raxsup=f"results/raxml/{output_prefix}/rxmlcatg.txt.raxml.support",
        # Use the same bamlist as the ANGSD/RAxML analysis (may include outgroups)
        bamlist=rules.make_bamlist_raxml_analysis.output.bamlist,
        samples=config["samples"],
        geno=rules.angsd_raxml.output.geno
    output:
        plots=f"{RAXML_FIG_DIR}/raxml_bootstrap.pdf"
    params:
        group_col = config["group_col"],
        populations = config["populations"]
    conda:
        "../envs/plot_tree.yaml"
    script:
        "../scripts/plot_raxml.R"


RAXML_TARGETS = [
    rules.angsd_raxml.output.geno,
    rules.angsd_raxml.output.bcf,
    f"results/angsd_raxml/{output_prefix}/gl.vcf.gz",
    rules.catg_format.output.raxmlcatg,
    rules.raxml_ng.output.raxout,
    rules.raxml_ng.output.raxsup,
    rules.plot_raxml.output.plots,
] if RAXML_ENABLED else []
