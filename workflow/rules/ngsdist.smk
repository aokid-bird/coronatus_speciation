# rules/ngsdist.smk

NGSDIST_CFG = (config.get("ngsdist", {}) or {})
_MODELS = NGSDIST_CFG.get("models", ["p", "jc69"])  # allowed: p, jc69
NGSDIST_SIF = f"{config_singularity_dir}/ngsdist.sif"

rule ngsdist_prepare_inputs:
    """
    Prepare labels and position info for ngsDist from ANGSD global outputs.
    - labels: one sample id per line (derived from BAM basenames)
    - posinfo: TSV with chr and site derived from beagle markers
    - nsites: number of variant sites (lines in beagle minus header)
    """
    input:
        bamlist = rules.angsd_global.input.bamlist,
        beagle  = rules.angsd_global.output.beagle
    output:
        labels = f"results/ngsdist_global/{output_prefix}/labels.txt",
        posinfo= f"results/ngsdist_global/{output_prefix}/posinfo.tsv",
        nsites = f"results/ngsdist_global/{output_prefix}/nsites.txt"
    conda:
        "../envs/intersect_sites.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.labels})
        # labels from bamlist -> sample ids
        awk -F'/' '{{print $NF}}' {input.bamlist} | sed -E 's/\.bam$//' > {output.labels}

        # posinfo from beagle markers (col1 like scaf_pos)
        python workflow/scripts/make_ngsld_snppos.py --beagle {input.beagle} --out {output.posinfo}

        # number of sites (lines - header)
        zcat {input.beagle} | tail -n +2 | wc -l | awk '{{print $1}}' > {output.nsites}
        """

rule ngsdist_run:
    """
    Run ngsDist with a selected evolutionary model: 'p' (0) or 'jc69' (1).
    Outputs are organized under .../ngsdist_global/{output_prefix}/{model}/.
    """
    input:
        beagle  = rules.angsd_global.output.beagle,
        labels  = rules.ngsdist_prepare_inputs.output.labels,
        posinfo = rules.ngsdist_prepare_inputs.output.posinfo,
        nsites  = rules.ngsdist_prepare_inputs.output.nsites,
        bamlist = rules.angsd_global.input.bamlist
    output:
        dist    = f"results/ngsdist_global/{output_prefix}/{{model}}/ngsdist",
        log     = f"logs/{output_prefix}/ngsdist_{{model}}.log"
    singularity:
        NGSDIST_SIF
    wildcard_constraints:
        model="|".join(_MODELS)
    shell:
        r"""
        # Map model name to ngsDist evol_model code
        EVOL=0
        if [ "{wildcards.model}" = "jc69" ]; then EVOL=1; fi

        # Use --pos only for JC69
        POS_OPT=""
        if [ "$EVOL" = "1" ]; then POS_OPT="--pos {input.posinfo}"; fi

        ngsDist \
          --probs true \
          --avg_nuc_dist \
          --evol_model $EVOL \
          --labels {input.labels} \
          $POS_OPT \
          --geno {input.beagle} \
          --n_ind $(wc -l < {input.bamlist}) \
          --n_sites $(cat {input.nsites}) \
          --out {output.dist} 2> {output.log}
        """

rule ngsdist_nexus:
    """
    Convert ngsDist distance matrix to NEXUS using model-specific directory.
    """
    input:
        dist = rules.ngsdist_run.output.dist
    output:
        nexus = f"results/ngsdist_global/{output_prefix}/{{model}}/ngsdist_input.nexus"
    conda:
        "../envs/ngsdist_post.yaml"
    wildcard_constraints:
        model="|".join(_MODELS)
    shell:
        """
        Rscript workflow/scripts/ngsdist_to_nexus.R {input.dist} {output.nexus} yes
        """
