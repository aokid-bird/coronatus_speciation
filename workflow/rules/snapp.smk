"""
SNAPP preparation workflow, including a dedicated ANGSD call and XML setup.
"""

SNAPP_MAX_PER_POP = int(SNAPP_CFG.get("max_per_population", 4))
SNAPP_MIN_SAMPLES_LOCUS = int(SNAPP_CFG.get("min_samples_locus", 4))
SNAPP_MISSINGNESS_THRESHOLD = (
    float(SNAPP_CFG.get("missingness_threshold"))
    if SNAPP_CFG.get("missingness_threshold") is not None
    else None
)
SNAPP_EXCLUDE_SAMPLES = _parse_list(SNAPP_CFG.get("exclude_samples"))
SNAPP_POPULATION_ALIASES = {
    str(key): str(value)
    for key, value in (SNAPP_CFG.get("population_aliases", {}) or {}).items()
}
SNAPP_OUTGROUP_ALIASES = {
    str(key): str(value)
    for key, value in (SNAPP_CFG.get("outgroup_aliases", {}) or {}).items()
}
SNAPP_CONSTRAINTS_CFG = SNAPP_CFG.get("constraints", {}) or {}
SNAPP_CONSTRAINT_TYPE = SNAPP_CONSTRAINTS_CFG.get("type") or SNAPP_CONSTRAINTS_CFG.get("placement")
_SNAPP_TAXA_CFG = SNAPP_CONSTRAINTS_CFG.get("taxa")
SNAPP_CONSTRAINT_TAXA = (
    ",".join(str(value) for value in _SNAPP_TAXA_CFG)
    if isinstance(_SNAPP_TAXA_CFG, (list, tuple))
    else _SNAPP_TAXA_CFG
)
SNAPP_RUN_MAP = SNAPP_CONSTRAINTS_CFG.get("runs", {}) or {}
SNAPP_RUN_IDS = sorted(SNAPP_RUN_MAP.keys())
SNAPP_PREP_CFG = SNAPP_CFG.get("snapp_prep", {}) or {}
SNAPP_MCMC_LENGTH = int(SNAPP_PREP_CFG.get("mcmc_length", 500000))
SNAPP_TOPOLOGY_WEIGHT = float(SNAPP_PREP_CFG.get("topology_weight", 1.0))
SNAPP_PREP_EXTRA = SNAPP_PREP_CFG.get("extra_args", "").strip()
SNAPP_LOG_PREFIX = SNAPP_CFG.get("log_prefix", "snapp")
SNAPP_ANGSD_MEM_MB = _resolve_mem_mb(16000, "analyses", "snapp", legacy_section=SNAPP_CFG)
SNAPP_ANGSD_RUNTIME = _resolve_runtime(1440, "analyses", "snapp", legacy_section=SNAPP_CFG)

SNAPP_BASE = f"results/snapp/{output_prefix}"
SNAPP_ALIGN_DIR = f"{SNAPP_BASE}/alignment"
SNAPP_VCF_DIR = f"{SNAPP_BASE}/vcf"
SNAPP_RUN_DIR = f"{SNAPP_BASE}/runs"
SNAPP_CONSTRAINT_DIR = f"{SNAPP_BASE}/constraints"


rule snapp_select_samples:
    """Select low-missingness individuals for SNAPP."""
    input:
        geno=f"results/angsd_global/{output_prefix}/gl.geno.gz",
        bamlist=f"results/bamlists/{output_prefix}/global_analysis/bamlist.txt",
        samples=config["samples"]
    output:
        summary=f"{SNAPP_BASE}/selected_samples.tsv",
        selected=f"{SNAPP_BASE}/selected_samples.txt"
    params:
        group_col=group_col,
        populations=" ".join(groups),
        max_per=SNAPP_MAX_PER_POP,
        optional_args=(
            " "
            + " ".join(
                (
                    ([f"--missingness-threshold {SNAPP_MISSINGNESS_THRESHOLD}"] if SNAPP_MISSINGNESS_THRESHOLD is not None else [])
                    + [f"--exclude-sample {s}" for s in SNAPP_EXCLUDE_SAMPLES]
                )
            )
            if (SNAPP_MISSINGNESS_THRESHOLD is not None or SNAPP_EXCLUDE_SAMPLES)
            else ""
        )
    conda:
        "../envs/python_pandas.yaml"
    shell:
        """
        python workflow/scripts/snapp_select_samples.py \
            --geno {input.geno} \
            --bamlist {input.bamlist} \
            --samples {input.samples} \
            --group-col {params.group_col} \
            --populations {params.populations} \
            --max-per-pop {params.max_per}{params.optional_args} \
            --summary-tsv {output.summary} \
            --selected-list {output.selected}
        """


rule snapp_make_bamlist:
    """Generate bamlist for SNAPP run including optional outgroups."""
    input:
        selected=rules.snapp_select_samples.output.selected,
        outgroups=(lambda wc: sliced_outgroup_inputs(SNAPP_OUTGROUP_IDS))
    output:
        bamlist=f"results/bamlists/{output_prefix}/snapp/bamlist.txt"
    params:
        bam_dir=config_bam_dir
    run:
        from pathlib import Path
        sel = Path(input.selected)
        samples = [line.strip() for line in sel.read_text().splitlines() if line.strip()]
        bam_paths = ingroup_bam_paths(sample_ids=samples, bam_dir=params.bam_dir)
        write_bamlist(output.bamlist, bam_paths + outgroup_bam_paths(SNAPP_OUTGROUP_IDS))


rule angsd_snapp:
    """ANGSD call restricted to SNAPP individuals."""
    input:
        bamlist=rules.snapp_make_bamlist.output.bamlist,
        sites=f"results/intersect_sites/{output_prefix}/intersect.txt",
        scafs=f"results/intersect_sites/{output_prefix}/intersect.chr",
        sites_idx=f"results/intersect_sites/{output_prefix}/intersect.txt.bin"
    output:
        geno=f"results/angsd_snapp/{output_prefix}/gl.geno.gz",
        mafs=f"results/angsd_snapp/{output_prefix}/gl.mafs.gz",
        beagle=f"results/angsd_snapp/{output_prefix}/gl.beagle.gz",
        bcf=f"results/angsd_snapp/{output_prefix}/gl.bcf"
    log:
        f"logs/{output_prefix}/angsd_snapp.log"
    params:
        ref=REF,
        outprefix=f"results/angsd_snapp/{output_prefix}/gl",
        extra=(config["angsd_common_args"].strip() + " " + (config.get("angsd_args", {}).get("snapp", "").strip())),
        minInd_ratio=get_minInd_ratio("snapp", get_minInd_ratio("global", None))
    threads: ANGSD_SNAPP_THREADS
    resources:
        mem_mb=SNAPP_ANGSD_MEM_MB,
        runtime=SNAPP_ANGSD_RUNTIME
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
              -sites {input.sites} -rf {input.scafs} \
              {params.extra} $MININD_OPT \
              -nThreads {threads} \
              2> {log}
        """


rule snapp_extract_unlinked:
    """
    Extract the LD-pruned site set used to restrict the SNAPP alignment.
    """
    input:
        beagle=rules.angsd_snapp.output.beagle,
        unlinked=f"results/ngsld_global/{output_prefix}/LD_unlinked.id"
    output:
        sites=f"{SNAPP_BASE}/unlinked_sites.tsv",
        summary=f"{SNAPP_BASE}/unlinked_summary.tsv"
    conda:
        "../envs/python_pandas.yaml"
    shell:
        """
        python workflow/scripts/snapp_extract_unlinked_sites.py \
            --beagle {input.beagle} \
            --unlinked-id {input.unlinked} \
            --out-sites {output.sites} \
            --out-summary {output.summary}
        """


rule snapp_prepare_metadata:
    """
    Build the sample-name and population metadata tables required by SNAPP helpers.
    """
    input:
        bamlist=rules.snapp_make_bamlist.output.bamlist,
        summary=rules.snapp_select_samples.output.summary,
        outgroups=config["outgroups"]
    output:
        sample_names=f"{SNAPP_BASE}/sample_names.txt",
        pop_table=f"{SNAPP_BASE}/populations.tsv",
        metadata=f"{SNAPP_BASE}/sample_metadata.tsv"
    params:
        group_col=group_col,
        pop_alias_args=" ".join([f"--population-alias {k}={v}" for k, v in SNAPP_POPULATION_ALIASES.items()]),
        out_alias_args=" ".join([f"--outgroup-alias {k}={v}" for k, v in SNAPP_OUTGROUP_ALIASES.items()])
    conda:
        "../envs/python_pandas.yaml"
    shell:
        """
        python workflow/scripts/snapp_prepare_metadata.py \
            --bamlist {input.bamlist} \
            --summary {input.summary} \
            --outgroups {input.outgroups} \
            --group-col {params.group_col} \
            --sample-names-out {output.sample_names} \
            --pop-table {output.pop_table} \
            --metadata-out {output.metadata} \
            {params.pop_alias_args} {params.out_alias_args}
        """


rule snapp_filter_vcf:
    """
    Filter the SNAPP BCF to the selected unlinked sites and reheader sample names.
    """
    input:
        bcf=rules.angsd_snapp.output.bcf,
        sites=rules.snapp_extract_unlinked.output.sites,
        sample_names=rules.snapp_prepare_metadata.output.sample_names
    output:
        vcf=f"{SNAPP_VCF_DIR}/snapp_unlinked.vcf.gz",
        tbi=f"{SNAPP_VCF_DIR}/snapp_unlinked.vcf.gz.tbi"
    log:
        f"logs/{output_prefix}/snapp_filter_vcf.log"
    conda:
        "../envs/bcftools_env.yaml"
    shell:
        """
        set -euo pipefail
        VCF_DIR=$(dirname {output.vcf})
        mkdir -p "$VCF_DIR"
        TMP_BCF=$(mktemp "$VCF_DIR/snapp_filter.XXXXXX.bcf")
        TMP_BCF_REHEADER="$TMP_BCF.rehead.bcf"
        cleanup() {{
            rm -f "$TMP_BCF" "$TMP_BCF_REHEADER"
        }}
        trap cleanup EXIT

        bcftools +fill-tags {input.bcf} -Ou -- -t HWE \
            | bcftools view -Ob -e 'HWE<=0.05' -T {input.sites} -o "$TMP_BCF"
        bcftools reheader -s {input.sample_names} -o "$TMP_BCF_REHEADER" "$TMP_BCF"
        bcftools view -Oz -o {output.vcf} "$TMP_BCF_REHEADER"
        bcftools index -t -f {output.vcf}
        """


rule snapp_vcf_to_phylip:
    """
    Convert the filtered SNAPP VCF into PHYLIP and binary NEXUS alignments.
    """
    input:
        vcf=rules.snapp_filter_vcf.output.vcf
    output:
        phylip=f"{SNAPP_ALIGN_DIR}/snapp.min{SNAPP_MIN_SAMPLES_LOCUS}.phy",
        nexus=f"{SNAPP_ALIGN_DIR}/snapp.min{SNAPP_MIN_SAMPLES_LOCUS}.bin.nexus"
    params:
        folder=SNAPP_ALIGN_DIR,
        min_samples=SNAPP_MIN_SAMPLES_LOCUS
    conda:
        "../envs/python_pandas.yaml"
    shell:
        """
        python workflow/scripts/vcf2phylip.py \
            --input {input.vcf} \
            --output-folder {params.folder} \
            --output-prefix snapp \
            --min-samples-locus {params.min_samples} \
            --nexus-binary
        """


rule snapp_constraints:
    """
    Write the per-run SNAPP topology constraint file from config settings.
    """
    output:
        constraints=f"{SNAPP_CONSTRAINT_DIR}/{{run}}.txt"
    params:
        type=SNAPP_CONSTRAINT_TYPE,
        taxa=SNAPP_CONSTRAINT_TAXA,
        dist_args=lambda wc: " ".join(f"--distribution '{d}'" for d in (SNAPP_RUN_MAP.get(wc.run, {}).get("distributions", []) or [])),
        dist_count=lambda wc: len(SNAPP_RUN_MAP.get(wc.run, {}).get("distributions", []) or [])
    conda:
        "../envs/python_pandas.yaml"
    shell:
        """
        if [ -z "{params.type}" ] || [ -z "{params.taxa}" ]; then
            echo "SNAPP constraints require 'type' and 'taxa' definitions" >&2
            exit 1
        fi
        if [ {params.dist_count} -eq 0 ]; then
            echo "No distributions defined for run {wildcards.run}" >&2
            exit 1
        fi
        python workflow/scripts/snapp_write_constraints.py \
            --type {params.type} \
            --taxa {params.taxa} \
            {params.dist_args} \
            --output {output.constraints}
        """


rule snapp_prep_xml:
    """
    Prepare the BEAST/SNAPP XML input for a configured constrained run.
    """
    input:
        phylip=rules.snapp_vcf_to_phylip.output.phylip,
        pop_table=rules.snapp_prepare_metadata.output.pop_table,
        constraints=f"{SNAPP_CONSTRAINT_DIR}/{{run}}.txt"
    output:
        xml=f"{SNAPP_RUN_DIR}/{{run}}/snapp.xml"
    params:
        length=SNAPP_MCMC_LENGTH,
        weight=SNAPP_TOPOLOGY_WEIGHT,
        extra=SNAPP_PREP_EXTRA,
        outprefix=f"{SNAPP_LOG_PREFIX}_{{run}}"
    log:
        f"logs/{output_prefix}/snapp_prep_{{run}}.log"
    conda:
        "../envs/snapp_prep.yaml"
    shell:
        """
        mkdir -p $(dirname {output.xml})
        ruby workflow/scripts/snapp_prep.rb \
            -p {input.phylip} \
            -t {input.pop_table} \
            -c {input.constraints} \
            -l {params.length} \
            -w {params.weight} \
            -o {params.outprefix} \
            -x {output.xml} \
            {params.extra} \
            > {log} 2>&1
        """


SNAPP_TARGETS = []
if SNAPP_ENABLED:
    SNAPP_TARGETS = [
        f"{SNAPP_BASE}/selected_samples.tsv",
        f"{SNAPP_BASE}/selected_samples.txt",
        f"results/bamlists/{output_prefix}/snapp/bamlist.txt",
        f"{SNAPP_BASE}/unlinked_sites.tsv",
        f"{SNAPP_BASE}/unlinked_summary.tsv",
        f"{SNAPP_BASE}/sample_names.txt",
        f"{SNAPP_BASE}/populations.tsv",
        f"{SNAPP_BASE}/sample_metadata.tsv",
        f"{SNAPP_VCF_DIR}/snapp_unlinked.vcf.gz",
        f"{SNAPP_VCF_DIR}/snapp_unlinked.vcf.gz.tbi",
        f"{SNAPP_ALIGN_DIR}/snapp.min{SNAPP_MIN_SAMPLES_LOCUS}.phy",
        f"{SNAPP_ALIGN_DIR}/snapp.min{SNAPP_MIN_SAMPLES_LOCUS}.bin.nexus",
    ]
    for run_id in SNAPP_RUN_IDS:
        SNAPP_TARGETS.extend([
            f"{SNAPP_BASE}/constraints/{run_id}.txt",
            f"{SNAPP_RUN_DIR}/{run_id}/snapp.xml",
        ])
