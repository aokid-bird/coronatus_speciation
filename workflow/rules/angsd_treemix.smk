###
# Treemix flow: ANGSD -> VCF -> ACF -> rename/meld/root -> treemix matrix -> block size -> treemix runs -> eval
###

GLACTOOLS_BIN = TREEMIX_CFG.get("glactools_bin", "workflow/bin/glactools")
TREEMIX_MODE = str(TREEMIX_CFG.get("mode", "merge"))
TREEMIX_ROOT_LABEL = str(TREEMIX_CFG.get("root_label", "")).strip() or None
TREEMIX_EXCLUDE_SAMPLES = _parse_list(TREEMIX_CFG.get("exclude_samples"))
TREEMIX_THREADS = _resolve_threads(TREEMIX_CFG, LEGACY_GLOBAL_THREADS)
TREEMIX_MAX_M = int(TREEMIX_CFG.get("max_m", 6))
TREEMIX_REPS = int(TREEMIX_CFG.get("reps", 10))
TREEMIX_TIMEOUT_SECONDS = int(TREEMIX_CFG.get("timeout_seconds", 300))
TREEMIX_MAX_ATTEMPTS = int(TREEMIX_CFG.get("max_attempts", 3))
_PLOTFUNC_DEFAULT = Path("workflow/scripts/treemix_plotting_funcs.R")
TREEMIX_PLOTTING_FUNCS = str(_PLOTFUNC_DEFAULT) if _PLOTFUNC_DEFAULT.exists() else None
TREEMIX_DIR = f"results/treemix/{output_prefix}/{TREEMIX_MODE}"
TREEMIX_PLOT_DIR = f"figures/exploratory/treemix/{output_prefix}/{TREEMIX_MODE}"

rule make_bamlist_treemix:
    """
    Write the TreeMix bamlist after applying per-analysis exclusions and outgroup selection.
    """
    input:
        ingroup=rules.make_bamlist_unrelated.output.bamlist,
        # Ensure sliced outgroup BAMs exist when including outgroups
        sliced=(lambda wc: sliced_outgroup_inputs(TREEMIX_OUTGROUP_IDS))
    output:
        bamlist=f"results/bamlists/{output_prefix}/treemix/bamlist.txt"
    run:
        import pandas as pd
        # start with ingroup bamlist
        ing = pd.read_csv(input.ingroup, header=None)[0].tolist()
        ing = filter_bam_paths_by_sample_ids(ing, exclude_samples=TREEMIX_EXCLUDE_SAMPLES)
        write_bamlist(output.bamlist, ing + outgroup_bam_paths(TREEMIX_OUTGROUP_IDS))


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
        acf=f"{TREEMIX_DIR}/input.acf.gz",
        treemix=f"{TREEMIX_DIR}/input.treemix"
    log:
        f"logs/{output_prefix}/glactools_prepare_treemix.log"
    params:
        group_col=group_col,
        glactools=GLACTOOLS_BIN,
        out_opt=(lambda wc: f"--include-outgroups --outgroups-tsv {config['outgroups']}" if TREEMIX_OUTGROUP_IDS else ""),
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
        f"{TREEMIX_DIR}/input.treemix"
    output:
        f"{TREEMIX_DIR}/input.treemix.gz"
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
        block=f"{TREEMIX_DIR}/lddecay/block_size.txt",
        ldplot=f"{TREEMIX_DIR}/lddecay/ld_plot.png"
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
        tm=f"{TREEMIX_DIR}/input.treemix.gz",
        block=f"{TREEMIX_DIR}/lddecay/block_size.txt"
    output:
        done=f"{TREEMIX_DIR}/edge_{{edge}}.done"
    params:
        # Build comma-delimited -root labels based on config/outgroups
        root_opt=(
            lambda wc: treemix_root_option(
                TREEMIX_OUTGROUP_IDS,
                merge_outgroups=TREEMIX_CFG.get("merge_outgroups", False),
                root_label=TREEMIX_ROOT_LABEL,
            )
        ),
        se_flag=(lambda wc: "-se" if (TREEMIX_CFG.get("se", True)) else ""),
        timeout_seconds=TREEMIX_TIMEOUT_SECONDS,
        max_attempts=TREEMIX_MAX_ATTEMPTS,
        reps=TREEMIX_REPS
    threads: 1
    conda:
        "../envs/treemix.yaml"
    shell:
        r"""
        EDGE={wildcards.edge}
        mkdir -p $(dirname {output.done})

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

        for REP in $(seq 1 {params.reps}); do
          OUTPREF={TREEMIX_DIR}/treemix_e${{EDGE}}_o${{REP}}
          LOG="${{OUTPREF}}.log"
          python workflow/scripts/run_treemix_with_timeout.py \
            --input {input.tm} \
            --outprefix "$OUTPREF" \
            --edge "$EDGE" \
            --block "$BLOCK" \
            --timeout-seconds {params.timeout_seconds} \
            --max-attempts {params.max_attempts} \
            --root-opt "{params.root_opt}" \
            --se-flag "{params.se_flag}" \
            --bootstrap \
            --log "$LOG"
        done

        touch {output.done}
        """


rule treemix_eval:
    input:
        expand(
            f"{TREEMIX_DIR}/edge_{{edge}}.done",
            edge=range(TREEMIX_MAX_M + 1),
        )
    output:
        summary=f"{TREEMIX_DIR}/eval_summary.csv",
        runs=f"{TREEMIX_DIR}/eval_runs.csv"
    params:
        runs_dir=TREEMIX_DIR,
        max_m=TREEMIX_MAX_M,
        out_prefix=f"{TREEMIX_DIR}/eval",
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
        runs=f"{TREEMIX_DIR}/eval_runs.csv",
        summary=f"{TREEMIX_DIR}/eval_summary.csv",
        bamlist=rules.make_bamlist_treemix.output.bamlist,
        samples=config["samples"]
    output:
        llkvar_pdf=f"{TREEMIX_PLOT_DIR}/treemix_validation2.pdf",
        resid_pdf=f"{TREEMIX_PLOT_DIR}/treemix_resid_{TREEMIX_MODE}.pdf",
        trees_pdf=f"{TREEMIX_PLOT_DIR}/treemix_trees_{TREEMIX_MODE}.pdf"
    params:
        runs_dir=TREEMIX_DIR,
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
        mkdir -p {TREEMIX_PLOT_DIR}
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


TREEMIX_TARGETS = [
    rules.angsd_treemix.output.geno,
    rules.angsd_treemix.output.bcf,
    rules.bcf_to_acf.output.vcf,
    rules.bcf_to_acf.output.acf,
    rules.glactools_prepare_treemix.output.acf,
    rules.glactools_prepare_treemix.output.treemix,
    f"{TREEMIX_DIR}/input.treemix.gz",
    rules.lddecay_blocksize.output.block,
    rules.lddecay_blocksize.output.ldplot,
    rules.treemix_eval.output.summary,
    rules.treemix_eval.output.runs,
    rules.treemix_plots.output.llkvar_pdf,
    rules.treemix_plots.output.resid_pdf,
    rules.treemix_plots.output.trees_pdf,
]
for _edge in range(TREEMIX_MAX_M + 1):
    TREEMIX_TARGETS.append(f"{TREEMIX_DIR}/edge_{_edge}.done")
if not TREEMIX_ENABLED:
    TREEMIX_TARGETS = []
