###
# Treemix flow: ANGSD -> VCF -> ACF -> rename/meld/root -> treemix matrix -> block size -> treemix runs -> eval
###

# Helper to construct a valid -root argument for TreeMix based on
# config outgroups and treemix settings, mirroring the RAxML-ng approach
# of deriving labels from available metadata.
def _abbr_taxon(t):
    try:
        import pandas as pd
        if pd.isna(t):
            return None
    except Exception:
        # if pandas not available in parser, fall back
        if t is None:
            return None
    words = str(t).strip().split()
    return "".join([w[:2] for w in words if len(w) > 0]) or None

def _treemix_root_opt(include_outgroups, mode, merge_outgroups, root_label):
    if not include_outgroups:
        return ""
    import pandas as pd
    og = pd.read_csv(config["outgroups"], sep="\t")
    # Optionally filter by species specified in config.treemix.outgroup_species
    species_cfg = (config.get("treemix", {}) or {}).get("outgroup_species", [])
    if isinstance(species_cfg, str):
        species_sel = [x.strip() for x in species_cfg.split(",") if x.strip()]
    else:
        species_sel = list(species_cfg or [])
    if species_sel and "taxon" in og.columns:
        og = og[og["taxon"].astype(str).isin(species_sel)].copy()

    labels = []
    if merge_outgroups:
        lab = root_label
        if not lab:
            if "taxon" in og.columns:
                ab = [x for x in (og["taxon"].apply(_abbr_taxon).dropna().tolist()) if x]
                uniq = sorted(list(set(ab)))
                lab = uniq[0] if len(uniq) == 1 else "OUT"
            else:
                lab = "OUT"
        labels = [lab]
    else:
        if "taxon" in og.columns:
            labels = sorted(list(set([x for x in og["taxon"].apply(_abbr_taxon).dropna().tolist() if x])))
        else:
            # fall back to per-sample labels
            labels = sorted(list(set(og["sample_id"].astype(str).tolist())))
    return "" if not labels else "-root " + ",".join(labels)

rule make_bamlist_treemix:
    """
    Bamlist for treemix, optionally including outgroups per config.treemix.include_outgroups.
    """
    input:
        ingroup=rules.make_bamlist_unrelated.output.bamlist,
        # Ensure sliced outgroup BAMs exist when including outgroups
        sliced=(lambda wc: [] if not TREEMIX_INCLUDE_OUTGROUPS else expand(f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam", sample_id=OUTGROUP_SAMPLE_IDS))
    output:
        bamlist=f"results/bamlists/{output_prefix}/treemix/bamlist.txt"
    params:
        bam_dir=config_bam_dir,
        include_out=TREEMIX_INCLUDE_OUTGROUPS
    run:
        import pandas as pd, os
        def _sample_id(path: str) -> str:
            name = os.path.basename(path)
            for suffix in (".bam", ".cram", ".sam"):
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
            if "_slice" in name:
                name = name.split("_slice", 1)[0]
            return name
        # start with ingroup bamlist
        ing = pd.read_csv(input.ingroup, header=None)[0].tolist()
        exclude = set(TREEMIX_EXCLUDE_SAMPLES)
        if exclude:
            ing = [path for path in ing if _sample_id(path) not in exclude]
        bams = list(ing)
        if params.include_out:
            for sid in OUTGROUP_SAMPLE_IDS:
                bams.append(f"{OUTGROUP_SLICED_DIR}/{sid}.bam")
            out_exclude = set(TREEMIX_EXCLUDE_OUTGROUPS)
            if out_exclude:
                og_prefix = f"{OUTGROUP_SLICED_DIR}/"
                bams = [
                    path
                    for path in bams
                    if not (path.startswith(og_prefix) and _sample_id(path) in out_exclude)
                ]
        pd.Series(bams).to_csv(output.bamlist, index=False, header=False)


rule angsd_treemix:
    """
    ANGSD on intersecting sites to produce BCF/geno for treemix.
    """
    input:
        bamlist=rules.make_bamlist_treemix.output.bamlist,
        sites=f"results/intersect_sites/{output_prefix}/intersect.txt",
        scafs=f"results/intersect_sites/{output_prefix}/intersect.chr",
        sites_idx = f"results/intersect_sites/{output_prefix}/intersect.txt.bin"
    output:
        geno=f"results/angsd_treemix/{output_prefix}/gl.geno.gz",
        bcf=f"results/angsd_treemix/{output_prefix}/gl.bcf"
    log:
        f"logs/{output_prefix}/angsd_treemix.log"
    params:
        ref=REF,
        outprefix=f"results/angsd_treemix/{output_prefix}/gl",
        extra=config["angsd_common_args"].strip() + " " + config["angsd_args"]["treemix"].strip(),
        minInd_ratio=get_minInd_ratio("treemix", get_minInd_ratio("global", None))
    threads: TREEMIX_THREADS
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


rule bcf_to_acf:
    """
    Convert ANGSD BCF to VCF and then to ACF using glactools binary.
    """
    input:
        bcf=f"results/angsd_treemix/{output_prefix}/gl.bcf",
        fai=FAI
    output:
        vcf=f"results/treemix/{output_prefix}/gl.vcf",
        acf=f"results/treemix/{output_prefix}/gl.acf.gz"
    log:
        f"logs/{output_prefix}/bcf_to_acf.log"
    params:
        prefix=f"results/treemix/{output_prefix}/gl",
        glactools=GLACTOOLS_BIN
    conda:
        "../envs/bcftools_env.yaml"
    shell:
        r"""
        bcftools convert -O v -o {params.prefix}.vcf {input.bcf}
        bcftools annotate -x FORMAT/GL {params.prefix}.vcf > {params.prefix}_simple.vcf
        {params.glactools} vcfm2acf --fai {input.fai} {params.prefix}_simple.vcf > {output.acf}
        {params.glactools} index {output.acf}
        if [ "{params.prefix}.vcf" != "{output.vcf}" ]; then mv {params.prefix}.vcf {output.vcf}; fi
        """


rule glactools_prepare_treemix:
    """
    Rename samples to IDs, optionally meld to populations, optionally set outgroup root, and write treemix matrix.
    """
    input:
        acf=rules.bcf_to_acf.output.acf,
        bamlist=rules.make_bamlist_treemix.output.bamlist,
        samples=config["samples"]
    output:
        acf=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/input.acf.gz",
        treemix=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/input.treemix"
    log:
        f"logs/{output_prefix}/glactools_prepare_treemix.log"
    params:
        group_col=group_col,
        glactools=GLACTOOLS_BIN,
        out_opt=(lambda wc: f"--include-outgroups --outgroups-tsv {config['outgroups']}" if TREEMIX_INCLUDE_OUTGROUPS else ""),
        merge_flag=(lambda wc: "true" if (TREEMIX_CFG.get("merge_outgroups", False)) else "false"),
        root_opt=(lambda wc: f"--root-label {TREEMIX_ROOT_LABEL}" if TREEMIX_ROOT_LABEL else ""),
        mode=TREEMIX_MODE
    conda:
        "../envs/python_pandas.yaml"
    shell:
        r"""
        python workflow/scripts/prepare_glactools_treemix.py \
            --acf {input.acf} \
            --bamlist {input.bamlist} \
            --samples {input.samples} \
            --group-col {params.group_col} \
            --mode {params.mode} \
            {params.out_opt} \
            --merge-outgroups {params.merge_flag} \
            {params.root_opt} \
            --glactools-bin {params.glactools} \
            --outdir $(dirname {output.acf}) \
            --treemix-input {output.treemix} \
            > {log} 2>&1
        """


rule treemix_input_gz:
    input:
        f"results/treemix/{output_prefix}/{TREEMIX_MODE}/input.treemix"
    output:
        f"results/treemix/{output_prefix}/{TREEMIX_MODE}/input.treemix.gz"
    shell:
        """
        gzip -c {input} > {output}
        """


rule lddecay_blocksize:
    input:
        ld=f"results/ngsld_global/{output_prefix}/LD.ld",
        acf=rules.glactools_prepare_treemix.output.acf,
        treemix=rules.glactools_prepare_treemix.output.treemix
    output:
        block=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/lddecay/block_size.txt",
        ldplot=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/lddecay/ld_plot.png"
    conda:
        "../envs/treemix_eval.yaml"
    shell:
        r"""
        export PATH=$(pwd)/workflow/bin:$PATH
        Rscript workflow/scripts/lddecay_blocksizefinder.R \
            --ld {input.ld} --acf {input.acf} --treemix {input.treemix} --outdir $(dirname {output.block})
        """


rule treemix_run:
    input:
        tm=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/input.treemix.gz",
        block=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/lddecay/block_size.txt"
    output:
        llik=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/treemix_e{{edge}}_o{{rep}}.llik"
    params:
        # Build comma-delimited -root labels based on config/outgroups
        root_opt=(
            lambda wc: _treemix_root_opt(
                TREEMIX_INCLUDE_OUTGROUPS,
                TREEMIX_MODE,
                TREEMIX_CFG.get("merge_outgroups", False),
                TREEMIX_ROOT_LABEL
            )
        ),
        se_flag=(lambda wc: "-se" if (TREEMIX_CFG.get("se", True)) else "")
    threads: 1
    conda:
        "../envs/treemix.yaml"
    shell:
        r"""
        EDGE={wildcards.edge}
        REP={wildcards.rep}
        OUTPREF=results/treemix/{output_prefix}/{TREEMIX_MODE}/treemix_e${{EDGE}}_o${{REP}}
        mkdir -p $(dirname {output.llik})

        # Validate block size; fallback to a sane default if needed
        BLOCK=$(cat {input.block} 2>/dev/null || echo "")
        if ! [[ "$BLOCK" =~ ^[0-9]+$ ]] || [ "$BLOCK" -le 0 ]; then
          echo "[treemix_run] Warning: invalid block size '$BLOCK'. Falling back to 1000." >&2
          BLOCK=1000
        fi
        # Ensure treemix is available
        if ! command -v treemix >/dev/null 2>&1; then
          echo "[treemix_run] Error: 'treemix' not found in PATH (check conda env)." >&2
          exit 127
        fi

        LOG="${{OUTPREF}}.log"
        # Pipe full output to .log without touching TreeMix's own .llik output file
        treemix -i {input.tm} -m $EDGE -o $OUTPREF {params.root_opt} -bootstrap -k $BLOCK {params.se_flag} 2>&1 | tee "$LOG"
        # Ensure the expected .llik file exists (written by TreeMix itself)
        test -s {output.llik}
        """


rule treemix_eval:
    input:
        expand(f"results/treemix/{output_prefix}/{TREEMIX_MODE}/treemix_e{{edge}}_o{{rep}}.llik", edge=range(TREEMIX_MAX_M+1), rep=range(1, TREEMIX_REPS+1))
    output:
        summary=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/eval_summary.csv",
        runs=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/eval_runs.csv"
    params:
        runs_dir=f"results/treemix/{output_prefix}/{TREEMIX_MODE}",
        max_m=TREEMIX_MAX_M,
        out_prefix=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/eval",
        plot_funcs_arg=(lambda wc: f"--plotting_funcs {TREEMIX_PLOTTING_FUNCS}" if TREEMIX_PLOTTING_FUNCS else ""),
        optm_flag=(lambda wc: "--optm" if TREEMIX_CFG.get("use_optm", True) else "")
    conda:
        "../envs/treemix_eval.yaml"
    shell:
        r"""
        Rscript workflow/scripts/treemix_optm_eval.R \
            --runs_dir {params.runs_dir} --max_m {params.max_m} --out_prefix {params.out_prefix} \
            {params.plot_funcs_arg} {params.optm_flag}
        # Move only if different to avoid 'same file' error
        if [ "{params.out_prefix}_summary.csv" != "{output.summary}" ]; then mv {params.out_prefix}_summary.csv {output.summary}; fi
        if [ "{params.out_prefix}_runs.csv"    != "{output.runs}" ];    then mv {params.out_prefix}_runs.csv {output.runs};    fi
        """


rule treemix_plots:
    """
    Generate TreeMix evaluation plots: lnL/VarExplain, residuals, and tree panels.
    """
    input:
        runs=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/eval_runs.csv",
        summary=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/eval_summary.csv",
        bamlist=rules.make_bamlist_treemix.output.bamlist,
        samples=config["samples"]
    output:
        llkvar_pdf=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/plots/treemix_validation2.pdf",
        resid_pdf=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/plots/treemix_resid_{TREEMIX_MODE}.pdf",
        trees_pdf=f"results/treemix/{output_prefix}/{TREEMIX_MODE}/plots/treemix_trees_{TREEMIX_MODE}.pdf"
    params:
        runs_dir=f"results/treemix/{output_prefix}/{TREEMIX_MODE}",
        mode=TREEMIX_MODE,
        group_col=group_col,
        plot_funcs=TREEMIX_PLOTTING_FUNCS,
        ggplot_treemix="workflow/scripts/ggplot_treemix.R",
        max_m=TREEMIX_MAX_M,
        reps=TREEMIX_REPS
    conda:
        "../envs/treemix_eval.yaml"
    shell:
        r"""
        unset R_LIBS_USER R_PROFILE_USER R_ENVIRON_USER
        Rscript workflow/scripts/treemix_make_plots.R \
            --runs_dir {params.runs_dir} \
            --mode {params.mode} \
            --samples {input.samples} \
            --group_col {params.group_col} \
            --bamlist {input.bamlist} \
            --max_m {params.max_m} \
            --reps {params.reps} \
            --outdir $(dirname {output.llkvar_pdf}) \
            --plotting_funcs {params.plot_funcs} \
            --ggplot_treemix {params.ggplot_treemix}
        """
