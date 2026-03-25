"""
Global ANGSD, relatedness, LD pruning, and unrelated-analysis inputs.
"""

kin_thr = config["ngsrelate"]["kinship_threshold"]
NGSRELATE_FIG_DIR = f"figures/exploratory/ngsrelate/{output_prefix}"
ANGSD_GLOBAL_MEM_MB = _resolve_mem_mb(16000, "analyses", "angsd_global", legacy_section=ANGSD_GLOBAL_CFG)
ANGSD_GLOBAL_RUNTIME = _resolve_runtime(1440, "analyses", "angsd_global", legacy_section=ANGSD_GLOBAL_CFG)
ANGSD_GLOBAL_UNRELATED_CFG = _config_section("angsd_global_unrelated_unlinked")
ANGSD_GLOBAL_UNRELATED_MEM_MB = _resolve_mem_mb(16000, "analyses", "angsd_global_unrelated_unlinked", legacy_section=ANGSD_GLOBAL_UNRELATED_CFG)
ANGSD_GLOBAL_UNRELATED_RUNTIME = _resolve_runtime(1440, "analyses", "angsd_global_unrelated_unlinked", legacy_section=ANGSD_GLOBAL_UNRELATED_CFG)

_INGROUP_BAMS = [
    f"{config_bam_dir}/{s}.bam"
    for s in ingroup_sample_ids(populations=groups)
]

rule make_bamlist_all:
    """
    Write the ingroup BAM list for the configured population set.
    """
    input:
        samples=config["samples"],
        bams=_INGROUP_BAMS
    output:
        bamlist=f"results/bamlists/{output_prefix}/global/bamlist.txt"
    run:
        write_bamlist(output.bamlist, ingroup_bam_paths(populations=groups))


rule make_bamlist_global_analysis:
    """
    Bamlist for angsd_global; optionally append sliced outgroups.
    """
    input:
        ingroup = rules.make_bamlist_all.output.bamlist,
        # Ensure sliced outgroup BAMs exist when including outgroups
        sliced=(lambda wc: sliced_outgroup_inputs(ANGSD_GLOBAL_OUTGROUP_IDS))
    output:
        bamlist = f"results/bamlists/{output_prefix}/global_analysis/bamlist.txt"
    run:
        import pandas as pd
        ing = pd.read_csv(input.ingroup, header=None)[0].tolist()
        write_bamlist(output.bamlist, list(ing) + outgroup_bam_paths(ANGSD_GLOBAL_OUTGROUP_IDS))

rule angsd_global:
    """
    Run ANGSD on the global analysis sample set and emit genotype-likelihood files.
    """
    input:
        bamlist=rules.make_bamlist_global_analysis.output.bamlist,
        sites=f"results/intersect_sites/{output_prefix}/intersect.txt",
        scafs=f"results/intersect_sites/{output_prefix}/intersect.chr",
        sites_idx = f"results/intersect_sites/{output_prefix}/intersect.txt.bin"
    output:
        geno=f"results/angsd_global/{output_prefix}/gl.geno.gz",
        mafs=f"results/angsd_global/{output_prefix}/gl.mafs.gz",
        beagle=f"results/angsd_global/{output_prefix}/gl.beagle.gz"
    log:
        f"logs/{output_prefix}/angsd_global.log"
    params:
        ref=REF,
        outprefix=f"results/angsd_global/{output_prefix}/gl",
        extra=config["angsd_common_args"].strip() + " " + config["angsd_args"]["global"].strip(),
        minInd_ratio=get_minInd_ratio("global", None)
    threads: ANGSD_GLOBAL_THREADS
    resources:
        mem_mb=ANGSD_GLOBAL_MEM_MB,
        runtime=ANGSD_GLOBAL_RUNTIME
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

rule make_ngsrelate_inputs:
    input:
        mafs   = f"results/angsd_global/{output_prefix}/gl.mafs.gz",
        bamlist=rules.make_bamlist_all.output.bamlist
    output:
        freq = f"results/ngsrelate_global/{output_prefix}/freq",
        id   = f"results/ngsrelate_global/{output_prefix}/id"
    params:
        bam_dir = config_bam_dir
    conda:
        "../envs/python_utils.yaml"
    shell:
        """
        python workflow/scripts/ngsrelate_prepare_inputs.py \
            --bamlist {input.bamlist} \
            --mafs {input.mafs} \
            --id_out {output.id} \
            --freq_out {output.freq} \
            --bam_dir {params.bam_dir}
        """

rule ngsrelate_global:
    input:
        beagle = f"results/angsd_global/{output_prefix}/gl.beagle.gz",
        freq   = rules.make_ngsrelate_inputs.output.freq,
        id     = rules.make_ngsrelate_inputs.output.id,
        bamlist= rules.make_bamlist_all.output.bamlist
    output:
        result = f"results/ngsrelate_global/{output_prefix}/ngsrelate_res"
    log:
        f"logs/{output_prefix}/ngsrelate_global.log"
    threads: NGSRELATE_THREADS
    singularity:
        f"{config_singularity_dir}/ngsrelate_20220925.sif"
    shell:
        """
        ngsRelate \
            -G {input.beagle} \
            -n $(wc -l < {input.bamlist}) \
            -f {input.freq} \
            -z {input.id} \
            -p {threads} \
            -O {output.result} > {log} 2>&1
        """

rule plot_ngsrelate_kinship:
    """
    Plot kinship summaries and write the related-sample table for manual filtering.
    """
    input:
        kin = f"results/ngsrelate_global/{output_prefix}/ngsrelate_res",
        meta = config["samples"]
    output:
        kinplot = f"{NGSRELATE_FIG_DIR}/kin_KING_thr{kin_thr}.pdf",
        netplot = f"{NGSRELATE_FIG_DIR}/kinnet_KING_thr{kin_thr}.pdf",
        csv = f"results/ngsrelate_global/{output_prefix}/related_KING_thr{kin_thr}.csv"
    params:
        kin_thr = config["ngsrelate"]["kinship_threshold"],
        pop_col = config["group_col"]
    conda:
        "../envs/r_plotting.yaml"
    shell:
        """
        mkdir -p {NGSRELATE_FIG_DIR}
        echo "R_LIBS_USER=${{R_LIBS_USER:-}}"
        echo "LD_LIBRARY_PATH=${{LD_LIBRARY_PATH:-}}"
        unset R_LIBS_USER R_PROFILE_USER R_ENVIRON_USER
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        Rscript workflow/scripts/plot_kinship.R {input.kin} {input.meta} {params.kin_thr} {params.pop_col}
        """

rule manual_remove_list:
    """
    Create the manual curation checkpoint for related-sample removal.
    """
    input:
        related_pairs = rules.plot_ngsrelate_kinship.output.csv
    output:
        remove_ids = f"results/ngsrelate_global/{output_prefix}/remove_ids.txt"
    params:
        kin_thr = config["ngsrelate"]["kinship_threshold"]
    message:
        "Please manually inspect related_KING_thr{params.kin_thr} and create {output.remove_ids}"

rule make_bamlist_unrelated:
    """
    Remove manually flagged related individuals from the global ingroup bamlist.
    """
    input:
        bamlist_global = rules.make_bamlist_all.output.bamlist,
        remove_ids = rules.manual_remove_list.output.remove_ids
    output:
        bamlist=f"results/bamlists/{output_prefix}/global_unrelated/bamlist.txt"
    run:
        import pandas as pd
        import os
        # load BAM list by lines
        bamlist = pd.read_csv(input.bamlist_global, header=None)[0]
        # load a list of ids to be removed
        try:
            remove_ids = pd.read_csv(input.remove_ids, header=None)[0].tolist()
        except pd.errors.EmptyDataError:# if no individuals to be removed
            remove_ids = []
        # Remove bam paths including remove ids
        bamlist_filtered = bamlist[~bamlist.apply(
            lambda path: os.path.splitext(os.path.basename(path))[0] in remove_ids
        )]
        write_bamlist(output.bamlist, bamlist_filtered.tolist())

rule make_ngsld_inputs:
    input:
        beagle = f"results/angsd_global/{output_prefix}/gl.beagle.gz"
    output:
        snppos = f"results/ngsld_global/{output_prefix}/snp.pos"
    conda:
        "../envs/python_utils.yaml"
    shell:
        """
        python workflow/scripts/make_ngsld_snppos.py \
            --beagle {input.beagle} \
            --out {output.snppos}
        """

rule ngsld_global:
    input:
        bamlist = rules.make_bamlist_global_analysis.output.bamlist,
        beagle  = f"results/angsd_global/{output_prefix}/gl.beagle.gz",
        snppos  = rules.make_ngsld_inputs.output.snppos
    output:
        ldout = f"results/ngsld_global/{output_prefix}/LD.ld",
    params:
        minmaf  = config["ngsld"]["minmaf"],
        maxdist = config["ngsld"]["max_kb_dist"],
        minmaf_flag = lambda wildcards: (
            f"--min_maf {config['ngsld']['minmaf']}"
            if config["ngsld"]["minmaf"] > 0 else ""
        )
    log:
        f"logs/{output_prefix}/ngsLD.log"
    threads:
        NGSLD_THREADS
    singularity:
        f"{config_singularity_dir}/ngsld_1.2.0.sif"
    shell:
        r"""
        ngsLD \
          --n_threads {threads} \
          --n_ind $(wc -l < {input.bamlist}) \
          --n_sites $(wc -l < {input.snppos}) \
          {params.minmaf_flag} \
          --ignore_miss_data \
          --probs --geno {input.beagle} \
          --pos   {input.snppos} \
          --max_kb_dist {params.maxdist} \
          --extend_out \
          --out {output.ldout} 2> {log}
        """

rule ld_pruning:
    input:
        ldout = rules.ngsld_global.output.ldout
    output:
        unlinkedid = f"results/ngsld_global/{output_prefix}/LD_unlinked.id"
    params:
        maxdist   = config["ngsld"]["max_kb_dist"] * 1000,
        minweight = config["ngsld"]["min_weight"],
        weight_field = config["ngsld"]["weight_field"]
    log:
        f"logs/{output_prefix}/pruneld.log"
    threads:
        NGSLD_THREADS
    singularity:
        f"{config_singularity_dir}/ngsld_1.2.0.sif"
    shell:
        r"""
        # remove NaNs
        awk 'BEGIN{{OFS="\t"}} !/NaN/' {input.ldout} > {input.ldout}.tmp
        mv {input.ldout}.tmp {input.ldout}

        # extract unlinked ids by prune_graph
        prune_graph \
          --n-threads {threads} \
          --in {input.ldout} \
          --weight-field "{params.weight_field}" \
          --weight-filter "column_3 <= {params.maxdist} && {params.weight_field} >= {params.minweight}" \
          --out {output.unlinkedid} 2> {log}
        """

rule filter_unlinked_and_summary:
    input:
        beagle    = f"results/angsd_global/{output_prefix}/gl.beagle.gz",
        unlinkedid= rules.ld_pruning.output.unlinkedid,
        snppos    = rules.make_ngsld_inputs.output.snppos
    output:
        unlinkedbeagle = f"results/ngsld_global/{output_prefix}/gl_global_unlinked.beagle.gz",
        summary        = f"results/ngsld_global/{output_prefix}/ld_summary.csv"
    conda:
        "../envs/python_utils.yaml"
    shell:
        """
        python workflow/scripts/process_unlinked_and_summary.py \
        --beagle {input.beagle} \
        --unlinkedid {input.unlinkedid} \
        --snppos {input.snppos} \
        --unlinkedbeagle {output.unlinkedbeagle} \
        --summary {output.summary}
        """

rule angsd_filtered_sites_index:
    """
    Convert LD-pruned site IDs into ANGSD-ready site and scaffold index files.
    """
    input:
        unlinkedid = rules.ld_pruning.output.unlinkedid
    output:
        sites      = f"results/unlinked_sites/{output_prefix}/sites.txt",
        scafs      = f"results/unlinked_sites/{output_prefix}/sites.chr",
        sites_idx  = f"results/unlinked_sites/{output_prefix}/sites.txt.bin"
    conda:
        "../envs/angsd.yaml"
    shell:
        r"""
        # 1) scaf:pos -> scaf\tpos, sort -> sites.txt
        awk -F':' '{{print $1"\t"$2}}' {input.unlinkedid} \
          | sort -k1,1 -k2,2n \
          > {output.sites}

        # 2) extract unique scafs from sites.txt -> sites.chr
        cut -f1 {output.sites} \
          | uniq \
          > {output.scafs}

        # 3) index sites -> sites.txt.bin
        angsd sites index {output.sites}
        """

rule make_bamlist_unrelated_analysis:
    """
    Write the unrelated-analysis bamlist, optionally appending sliced outgroups.
    """
    input:
        ingroup_unrel = rules.make_bamlist_unrelated.output.bamlist,
        # Ensure sliced outgroup BAMs exist when including outgroups
        sliced=(lambda wc: sliced_outgroup_inputs(ANGSD_GLOBAL_OUTGROUP_IDS))
    output:
        bamlist = f"results/bamlists/{output_prefix}/global_unrelated_analysis/bamlist.txt"
    run:
        import pandas as pd
        ing = pd.read_csv(input.ingroup_unrel, header=None)[0].tolist()
        write_bamlist(output.bamlist, list(ing) + outgroup_bam_paths(ANGSD_GLOBAL_OUTGROUP_IDS))


rule angsd_global_unrelated_unlinked:
    """
    Run ANGSD on unrelated individuals at the LD-pruned unlinked sites.
    """
    input:
        bamlist=rules.make_bamlist_unrelated_analysis.output.bamlist,
        sites=rules.angsd_filtered_sites_index.output.sites,
        scafs=rules.angsd_filtered_sites_index.output.scafs,
        sites_idx = rules.angsd_filtered_sites_index.output.sites_idx
    output:
        geno=f"results/angsd_global_unrelated_unlinked/{output_prefix}/gl.geno.gz",
        mafs=f"results/angsd_global_unrelated_unlinked/{output_prefix}/gl.mafs.gz",
        beagle=f"results/angsd_global_unrelated_unlinked/{output_prefix}/gl.beagle.gz"
    log:
        f"logs/{output_prefix}/angsd_global_unrelated_unlinked.log"
    params:
        ref=REF,
        outprefix=f"results/angsd_global_unrelated_unlinked/{output_prefix}/gl",
        extra=config["angsd_common_args"].strip() + " " + config["angsd_args"]["global_unrelated_unlinked"].strip(),
        minInd_ratio=get_minInd_ratio("global", None)
    threads: ANGSD_GLOBAL_UNRELATED_THREADS
    resources:
        mem_mb=ANGSD_GLOBAL_UNRELATED_MEM_MB,
        runtime=ANGSD_GLOBAL_UNRELATED_RUNTIME
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


ANGSD_GLOBAL_CORE_TARGETS = [
    f"results/angsd_global/{output_prefix}/gl.geno.gz" if ANGSD_GLOBAL_ACTIVE else [],
    f"results/ngsrelate_global/{output_prefix}/ngsrelate_res" if NGSRELATE_ACTIVE else [],
    f"results/ngsrelate_global/{output_prefix}/related_KING_thr{kin_thr}.csv" if NGSRELATE_ACTIVE else [],
    f"results/ngsld_global/{output_prefix}/snp.pos" if NGSLD_ACTIVE else [],
    f"results/ngsld_global/{output_prefix}/LD.ld" if NGSLD_ACTIVE else [],
    f"results/ngsld_global/{output_prefix}/LD_unlinked.id" if NGSLD_ACTIVE else [],
    f"results/ngsld_global/{output_prefix}/gl_global_unlinked.beagle.gz" if NGSLD_ACTIVE else [],
    [
        f"results/unlinked_sites/{output_prefix}/sites.txt",
        f"results/unlinked_sites/{output_prefix}/sites.chr",
        f"results/unlinked_sites/{output_prefix}/sites.txt.bin",
    ] if (SFS_ENABLED and SFS_NEEDS_UNLINKED_SITES) else [],
    f"results/angsd_global_unrelated_unlinked/{output_prefix}/gl.beagle.gz" if GLOBAL_UNRELATED_UNLINKED_ACTIVE else [],
]
