# rules/angsd_raxml.smk

# Helper: build --outgroup from full species names defined in config.
def _outgroup_opt_from_species(bamlist_path, species_list):
    import os
    import pandas as pd
    try:
        with open(bamlist_path) as f:
            names = [os.path.basename(l.strip()).replace(".bam", "") for l in f if l.strip()]
    except FileNotFoundError:
        return ""

    # Normalize species list: allow list or comma-delimited string
    if isinstance(species_list, str):
        wanted = [x.strip() for x in species_list.split(",") if x.strip()]
    else:
        wanted = list(species_list or [])

    # If user didn't specify, fall back to all outgroup sample IDs
    if not wanted:
        candidates = set(OUTGROUP_SAMPLE_IDS)
    else:
        # Read outgroups TSV and select sample_ids whose taxon matches the wanted species
        try:
            og = pd.read_csv(config["outgroups"], sep="\t")
            if "taxon" in og.columns:
                candidates = set(og.loc[og["taxon"].astype(str).isin(wanted), "sample_id"].astype(str).tolist())
            else:
                candidates = set()  # taxon column missing; no match possible
        except Exception:
            candidates = set()

    # Intersect with actual MSA labels (bam basenames)
    selected = [n for n in names if n in candidates]
    return "" if not selected else "--outgroup " + ",".join(selected)

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
        sliced=(lambda wc: [] if not ANGSD_RAXML_INCLUDE_OUTGROUPS else expand(f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam", sample_id=OUTGROUP_SAMPLE_IDS))
    output:
        bamlist = f"results/bamlists/{output_prefix}/raxml_analysis/bamlist.txt"
    params:
        include_out = ANGSD_RAXML_INCLUDE_OUTGROUPS
    run:
        import pandas as pd
        ing = pd.read_csv(input.ingroup_downsampled, header=None)[0].tolist()
        bams = list(ing)
        if params.include_out:
            bams += [f"{OUTGROUP_SLICED_DIR}/{sid}.bam" for sid in OUTGROUP_SAMPLE_IDS]
        pd.Series(bams).to_csv(output.bamlist, index=False, header=False)


rule angsd_raxml:
    """
    ANGSD using all samples of the selected populations 
    to a generate RAxML-ng input file while excluding too closely related individuals
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
    threads: config["threads"]
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
    Generating CATG input files from a vcf.gz output of ANGSD
    """
    input:
        vcf=f"results/angsd_raxml/{output_prefix}/gl.vcf.gz",
        bamlist=rules.make_bamlist_raxml_analysis.output.bamlist
    output:
        raxmlcatg=f"results/raxml/{output_prefix}/rxmlcatg.txt"
    params:
        maxSize = 4 * 1024**3 # max size for future.apply
    threads:config['raxml']['threads_catg']
    conda:
        "../envs/vcfR.yaml"
    script:
        "../scripts/catg_formatter.R"

rule raxml_ng:
    """
    Run RAxML-ng using CATG input format
    """
    input:
        raxmlcatg=rules.catg_format.output.raxmlcatg,
        # also pass bamlist to derive exact CATG labels for outgroup option
        bamlist=rules.make_bamlist_raxml_analysis.output.bamlist
    output:
        raxout=f"results/raxml/{output_prefix}/rxmlcatg.txt.raxml.bootstraps",
        raxsup=f"results/raxml/{output_prefix}/rxmlcatg.txt.raxml.support"
    threads:config['raxml']['threads_run']
    params:
        model=config['raxml']['model'],
        bs=config['raxml']['bs'],
        # Pre-render optional outgroup argument from actual bamlist labels
        outgroup_opt=(
            (lambda wildcards, input: _outgroup_opt_from_species(
                input.bamlist,
                (config.get('raxml', {}) or {}).get('outgroup_species', [])
             )) if ANGSD_RAXML_INCLUDE_OUTGROUPS else ""
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
        plots=f"figures/exploratory/raxml/{output_prefix}/raxml_bootstrap.pdf"
    params:
        group_col = config["group_col"],
        populations = config["populations"]
    conda:
        "../envs/plot_tree.yaml"
    script:
        "../scripts/plot_raxml.R"
