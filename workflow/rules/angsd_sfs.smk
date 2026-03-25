"""
Site-frequency spectrum, theta, Tajima's D, and pairwise Fst analyses.
"""

from pathlib import Path

ANGSD_SFS_MEM_MB = _resolve_mem_mb(16000, "analyses", "angsd_sfs", legacy_section=ANGSD_SFS_CFG)
ANGSD_SFS_RUNTIME = _resolve_runtime(1440, "analyses", "angsd_sfs", legacy_section=ANGSD_SFS_CFG)
REALSFS_CFG = _config_section("sfs_analysis", "realSFS")
REALSFS_MEM_MB = _resolve_mem_mb(16000, "analyses", "sfs_realsfs", legacy_section=REALSFS_CFG)
REALSFS_RUNTIME = _resolve_runtime(480, "analyses", "sfs_realsfs", legacy_section=REALSFS_CFG)

SFS_SITE_INPUTS = {
    "linked": {
        "sites": f"results/intersect_sites/{output_prefix}/intersect.txt",
        "scafs": f"results/intersect_sites/{output_prefix}/intersect.chr",
        "sites_idx": f"results/intersect_sites/{output_prefix}/intersect.txt.bin",
    },
    "unlinked": {
        "sites": f"results/unlinked_sites/{output_prefix}/sites.txt",
        "scafs": f"results/unlinked_sites/{output_prefix}/sites.chr",
        "sites_idx": f"results/unlinked_sites/{output_prefix}/sites.txt.bin",
    },
}

SUPPORTED_SITE_KEYS = set(SFS_SITE_INPUTS)
invalid_sites = set(SFS_SITE_FILTERS) - SUPPORTED_SITE_KEYS
if invalid_sites:
    raise ValueError(f"Invalid site filter(s) for SFS analysis: {sorted(invalid_sites)}")
SFS_SITE_FILTERS_ACTIVE = [sf for sf in SFS_SITE_FILTERS if sf in SUPPORTED_SITE_KEYS]
if not SFS_SITE_FILTERS_ACTIVE:
    SFS_SITE_FILTERS_ACTIVE = ["linked"]

SUPPORTED_FOLD_KEYS = {"fold", "unfold"}
invalid_folds = set(SFS_FOLD_STATES) - SUPPORTED_FOLD_KEYS
if invalid_folds:
    raise ValueError(f"Invalid fold state(s) for SFS analysis: {sorted(invalid_folds)}")
SFS_FOLD_STATES_ACTIVE = [fl for fl in SFS_FOLD_STATES if fl in SUPPORTED_FOLD_KEYS]
if not SFS_FOLD_STATES_ACTIVE:
    SFS_FOLD_STATES_ACTIVE = ["fold"]

PAIR_DELIM = "__"


def _pair_label(pop1, pop2):
    return f"{pop1}{PAIR_DELIM}{pop2}"


SFS_PAIR_LABELS = [_pair_label(a, b) for a, b in SFS_PAIRWISE_COMBOS]
SFS_PAIR_LOOKUP = {_pair_label(a, b): (a, b) for a, b in SFS_PAIRWISE_COMBOS}


def _ensure_pair(pair_label: str):
    try:
        return SFS_PAIR_LOOKUP[pair_label]
    except KeyError as exc:
        raise ValueError(f"Pair '{pair_label}' not configured for SFS 2D analysis") from exc


def _site_path(site_filter: str, key: str):
    try:
        return SFS_SITE_INPUTS[site_filter][key]
    except KeyError as exc:
        raise ValueError(f"Unsupported site filter '{site_filter}' for SFS analysis") from exc


def _angsd_sfs_prefix(site_filter: str, group: str) -> str:
    return f"results/angsd_sfs/{output_prefix}/{site_filter}/{group}/gl"


def _realsfs1d_prefix(site_filter: str, fold_state: str, group: str) -> str:
    return f"results/realsfs_1d/{output_prefix}/{site_filter}/{fold_state}/{group}/gl"


def _realsfs2d_prefix(site_filter: str, fold_state: str, pair: str) -> str:
    return f"results/realsfs_2d/{output_prefix}/{site_filter}/{fold_state}/{pair}/pair"


def _angsd_sfs_group_input_bamlist(wildcards):
    return f"results/bamlists/{output_prefix}/{wildcards.group}/bamlist.txt"


def _angsd_sfs_group_input_sites(wildcards):
    return _site_path(wildcards.site_filter, "sites")


def _angsd_sfs_group_input_scafs(wildcards):
    return _site_path(wildcards.site_filter, "scafs")


def _angsd_sfs_group_input_sites_idx(wildcards):
    return _site_path(wildcards.site_filter, "sites_idx")


def _angsd_sfs_group_param_outprefix(wildcards):
    return _angsd_sfs_prefix(wildcards.site_filter, wildcards.group)


def _realsfs1d_input_saf_idx(wildcards):
    return f"{_angsd_sfs_prefix(wildcards.site_filter, wildcards.group)}.saf.idx"


def _realsfs1d_param_outprefix(wildcards):
    return _realsfs1d_prefix(wildcards.site_filter, wildcards.fold, wildcards.group)


def _realsfs2d_input_saf1(wildcards):
    pair = _ensure_pair(wildcards.pair)
    return f"{_angsd_sfs_prefix(wildcards.site_filter, pair[0])}.saf.idx"


def _realsfs2d_input_saf2(wildcards):
    pair = _ensure_pair(wildcards.pair)
    return f"{_angsd_sfs_prefix(wildcards.site_filter, pair[1])}.saf.idx"


def _realsfs2d_param_outprefix(wildcards):
    return _realsfs2d_prefix(wildcards.site_filter, wildcards.fold, wildcards.pair)


rule angsd_sfs_group:
    """Run ANGSD per population and site filter to generate SAF files for SFS."""
    input:
        bamlist=_angsd_sfs_group_input_bamlist,
        sites=_angsd_sfs_group_input_sites,
        scafs=_angsd_sfs_group_input_scafs,
        sites_idx=_angsd_sfs_group_input_sites_idx,
    output:
        geno=f"results/angsd_sfs/{output_prefix}/{{site_filter}}/{{group}}/gl.geno.gz",
        saf_idx=f"results/angsd_sfs/{output_prefix}/{{site_filter}}/{{group}}/gl.saf.idx",
        saf=f"results/angsd_sfs/{output_prefix}/{{site_filter}}/{{group}}/gl.saf.gz",
        mafs=f"results/angsd_sfs/{output_prefix}/{{site_filter}}/{{group}}/gl.mafs.gz",
        bcf=f"results/angsd_sfs/{output_prefix}/{{site_filter}}/{{group}}/gl.bcf",
    log:
        f"logs/{output_prefix}/angsd_sfs_{{site_filter}}_{{group}}.log",
    params:
        outprefix=_angsd_sfs_group_param_outprefix,
        ref=REF,
        extra=config["angsd_common_args"].strip() + " " + config["angsd_args"]["sfs"].strip(),
        minInd_ratio=get_minInd_ratio("sfs", None),
    threads:
        ANGSD_SFS_THREADS
    resources:
        mem_mb=ANGSD_SFS_MEM_MB,
        runtime=ANGSD_SFS_RUNTIME
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {params.outprefix})

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


rule realsfs_1d:
    """Estimate 1D SFS, theta, and Tajima's D per population."""
    input:
        saf_idx=_realsfs1d_input_saf_idx,
    output:
        sfs=f"results/realsfs_1d/{output_prefix}/{{site_filter}}/{{fold}}/{{group}}/gl.sfs",
        thetas_idx=f"results/realsfs_1d/{output_prefix}/{{site_filter}}/{{fold}}/{{group}}/gl.thetas.idx",
        thetas=f"results/realsfs_1d/{output_prefix}/{{site_filter}}/{{fold}}/{{group}}/gl.thetas.gz",
        pest=f"results/realsfs_1d/{output_prefix}/{{site_filter}}/{{fold}}/{{group}}/gl.thetas.idx.pestPG",
    params:
        fold=lambda wc: 1 if wc.fold == "fold" else 0,
        outprefix=_realsfs1d_param_outprefix,
    log:
        f"logs/{output_prefix}/realsfs_1d_{{site_filter}}_{{fold}}_{{group}}.log",
    threads:
        REALSFS_THREADS
    resources:
        mem_mb=REALSFS_MEM_MB,
        runtime=REALSFS_RUNTIME
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.sfs})

        FOLD_OPT=""
        if [ {params.fold} -eq 1 ]; then
          FOLD_OPT="-fold 1"
        fi

        realSFS {input.saf_idx} -maxIter {SFS_MAXITER} -tole {SFS_TOLE} -P {threads} $FOLD_OPT > {output.sfs} 2> {log}
        realSFS saf2theta {input.saf_idx} -outname {params.outprefix} -sfs {output.sfs} $FOLD_OPT >> {log} 2>&1
        thetaStat do_stat {output.thetas_idx} >> {log} 2>&1
        """


rule realsfs_2d:
    """Compute 2D SFS and Fst across population pairs."""
    input:
        saf1=_realsfs2d_input_saf1,
        saf2=_realsfs2d_input_saf2,
    output:
        sfs=f"results/realsfs_2d/{output_prefix}/{{site_filter}}/{{fold}}/{{pair}}/pair.ml",
        fst_idx=f"results/realsfs_2d/{output_prefix}/{{site_filter}}/{{fold}}/{{pair}}/pair.fst.idx",
        fst_stats=f"results/realsfs_2d/{output_prefix}/{{site_filter}}/{{fold}}/{{pair}}/pair.fst.txt",
        fst_window=f"results/realsfs_2d/{output_prefix}/{{site_filter}}/{{fold}}/{{pair}}/pair.fst.window.tsv",
        bins=f"results/realsfs_2d/{output_prefix}/{{site_filter}}/{{fold}}/{{pair}}/pair.bin",
    params:
        fold=lambda wc: 1 if wc.fold == "fold" else 0,
        outprefix=_realsfs2d_param_outprefix,
    log:
        f"logs/{output_prefix}/realsfs_2d_{{site_filter}}_{{fold}}_{{pair}}.log",
    threads:
        REALSFS_THREADS
    resources:
        mem_mb=REALSFS_MEM_MB,
        runtime=REALSFS_RUNTIME
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.sfs})

        FOLD_OPT=""
        if [ {params.fold} -eq 1 ]; then
          FOLD_OPT="-fold 1"
        fi

        realSFS {input.saf1} {input.saf2} -maxIter {SFS_MAXITER} -tole {SFS_TOLE} -P {threads} $FOLD_OPT > {output.sfs} 2> {log}
        realSFS fst index {input.saf1} {input.saf2} -sfs {output.sfs} -fstout {params.outprefix} $FOLD_OPT >> {log} 2>&1
        realSFS fst stats {output.fst_idx} > {output.fst_stats} 2>> {log}
        realSFS bins {input.saf1} {input.saf2} -P {threads} $FOLD_OPT > {output.bins} 2>> {log}
        realSFS fst stats2 {output.fst_idx} -win 50000 -step 10000 > {output.fst_window} 2>> {log}
        """


rule summarize_sfs_stats:
    """Aggregate SFS-derived diversity statistics and run linear models."""
    input:
        sfs=[f"{_realsfs1d_prefix(site, fold, group)}.sfs" for site in SFS_SITE_FILTERS_ACTIVE for fold in SFS_FOLD_STATES_ACTIVE for group in groups],
        pest=[f"{_realsfs1d_prefix(site, fold, group)}.thetas.idx.pestPG" for site in SFS_SITE_FILTERS_ACTIVE for fold in SFS_FOLD_STATES_ACTIVE for group in groups],
        contigs=f"results/reference/{output_prefix}/contig_lengths.csv"
    output:
        sfs_summary=f"results/sfs/{output_prefix}/sfs_1d_summary.csv",
        theta_summary=f"results/sfs/{output_prefix}/theta_window_summary.csv",
        lm_summary=f"results/sfs/{output_prefix}/lm_theta_tajima_summary.csv",
        plot=f"results/sfs/{output_prefix}/tajima_vs_theta.pdf"
    params:
        include_monomorphic=SFS_INCLUDE_MONOMORPHIC,
        site_filters=SFS_SITE_FILTERS_ACTIVE,
        fold_states=SFS_FOLD_STATES_ACTIVE,
        populations=groups,
        min_contig_length=2_000_000
    conda:
        "../envs/r_plotting.yaml"
    script:
        "../scripts/summarize_sfs.R"


ANGSD_SFS_TARGETS = [
    *expand(
        f"results/angsd_sfs/{output_prefix}/{{site_filter}}/{{group}}/gl.geno.gz",
        site_filter=SFS_SITE_FILTERS_ACTIVE,
        group=groups,
    ),
    *expand(
        f"results/realsfs_1d/{output_prefix}/{{site_filter}}/{{fold}}/{{group}}/gl.sfs",
        site_filter=SFS_SITE_FILTERS_ACTIVE,
        fold=SFS_FOLD_STATES_ACTIVE,
        group=groups,
    ),
    *expand(
        f"results/realsfs_1d/{output_prefix}/{{site_filter}}/{{fold}}/{{group}}/gl.thetas.idx.pestPG",
        site_filter=SFS_SITE_FILTERS_ACTIVE,
        fold=SFS_FOLD_STATES_ACTIVE,
        group=groups,
    ),
    f"results/reference/{output_prefix}/contig_lengths.csv",
    *(
        expand(
            f"results/realsfs_2d/{output_prefix}/{{site_filter}}/{{fold}}/{{pair}}/pair.ml",
            site_filter=SFS_SITE_FILTERS_ACTIVE,
            fold=SFS_FOLD_STATES_ACTIVE,
            pair=SFS_PAIR_LABELS,
        )
        if SFS_PAIR_LABELS else []
    ),
    f"results/sfs/{output_prefix}/sfs_1d_summary.csv",
    f"results/sfs/{output_prefix}/theta_window_summary.csv",
    f"results/sfs/{output_prefix}/lm_theta_tajima_summary.csv",
    f"results/sfs/{output_prefix}/tajima_vs_theta.pdf",
] if SFS_ENABLED else []
