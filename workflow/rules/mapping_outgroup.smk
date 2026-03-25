"""
Mapping trimmed outgroup reads to the reference, with optional inclusion of
unpaired reads. Mapper can be selected via config.mapping.outgroup.mapper (or
fallback to the global mapping defaults): bwa or bwa-mem2.

Outputs final coordinate-sorted BAM and index in config['bam_dir']/outgroups/.
"""

import os
import re
from os.path import join as pjoin

OUTGROUP_MAP_TMP = OUTGROUP_MAP_TMP_DIR
OUTGROUP_FINAL_DIR = OUTGROUP_BAM_DIR

MAP_CFG = MAPPING_OUTGROUP_CFG
LR_CFG = (config.get("longread", {}) or {})
MAP_UNPAIRED = bool(MAP_CFG.get("map_unpaired", True))
MAP_EXTRA = str(MAP_CFG.get("extra", "") or "").strip()
OUTGROUP_MAP_PAIRED_MEM_MB = _resolve_mem_mb(64000, "mapping", "outgroup_paired")
OUTGROUP_MAP_PAIRED_RUNTIME = _resolve_runtime(1440, "mapping", "outgroup_paired")
OUTGROUP_MAP_UNPAIRED_MEM_MB = _resolve_mem_mb(64000, "mapping", "outgroup_unpaired")
OUTGROUP_MAP_UNPAIRED_RUNTIME = _resolve_runtime(1440, "mapping", "outgroup_unpaired")
OUTGROUP_MAP_FINAL_MEM_MB = _resolve_mem_mb(200000, "mapping", "outgroup_merge")
OUTGROUP_MAP_FINAL_RUNTIME = _resolve_runtime(1440, "mapping", "outgroup_merge")
OUTGROUP_LONGREAD_MAP_MEM_MB = _resolve_mem_mb(200000, "mapping", "outgroup_longread")
OUTGROUP_LONGREAD_MAP_RUNTIME = _resolve_runtime(1440, "mapping", "outgroup_longread")
OUTGROUP_FILTLONG_CFG = (LR_CFG.get("filter", {}) or {})
OUTGROUP_FILTLONG_THREADS = _resolve_threads(OUTGROUP_FILTLONG_CFG, 4, legacy_fallback=False)

# choose mapper subcommand and reference argument (prefix or fasta)
MAPPER_CMD = "bwa-mem2 mem" if MAPPER_OUTGROUP == "bwa-mem2" else "bwa mem"
MAPPER_REF_ARG = REF_MAP_ARG_OUTGROUP

# Derive sample ID subsets
LONGREAD_SAMPLES = [sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "long"]
SHORTREAD_SAMPLES = [sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "short"]

def _wc_regex(ids):
    # Build a regex that matches only the provided IDs; if empty, match nothing
    return "(" + "|".join(map(re.escape, ids)) + ")" if ids else r"a^"

# Optional: mock reference/index only for dry-run (no-op in real runs)
if config.get("dryrun_mock_reference", False):
    rule mock_reference_for_dryrun:
        """
        Expose reference fasta and mapper index as buildable outputs so dry-run can resolve the DAG
        even if these files are already present on disk in production.
        """
        output:
            ref=REF_UNZIPPED,
            mapidx=REF_MAP_INDEX_OUTGROUP,
            fai=FAI
        shell:
            """
            true
            """

rule map_outgroup_paired:
    """
    Map paired trimmed outgroup short reads with the configured short-read mapper.
    """
    wildcard_constraints:
        sample_id=_wc_regex(SHORTREAD_SAMPLES)
    input:
        ref=lambda wc: REF_UNZIPPED,
        idx=lambda wc: REF_MAP_INDEX_OUTGROUP,
        fq1=lambda wc: pjoin(OUTGROUP_TRIM_DIR, f"{wc.sample_id}_pair_R1.fastq.gz"),
        fq2=lambda wc: pjoin(OUTGROUP_TRIM_DIR, f"{wc.sample_id}_pair_R2.fastq.gz")
    output:
        bam=temp(pjoin(OUTGROUP_MAP_TMP, "{sample_id}.paired.bam"))
    threads: MAPPING_OUTGROUP_THREADS
    resources:
        mem_mb=OUTGROUP_MAP_PAIRED_MEM_MB,
        runtime=OUTGROUP_MAP_PAIRED_RUNTIME
    conda:
        "../envs/mapper.yaml"
    message:
        "Map paired reads for outgroup {wildcards.sample_id} using {MAPPER_CMD}"
    shell:
        r"""
        mkdir -p {OUTGROUP_MAP_TMP}
        {MAPPER_CMD} -t {threads} {MAP_EXTRA} {MAPPER_REF_ARG} {input.fq1} {input.fq2} \
            | samtools view -Sb - \
            > {output.bam}
        """

rule map_outgroup_unpaired:
    """
    Map trimmed outgroup unpaired reads when unpaired mapping is enabled.
    """
    wildcard_constraints:
        sample_id=_wc_regex(SHORTREAD_SAMPLES)
    input:
        ref=lambda wc: REF_UNZIPPED,
        idx=lambda wc: REF_MAP_INDEX_OUTGROUP,
        fq=lambda wc: pjoin(OUTGROUP_TRIM_DIR, f"{wc.sample_id}_unpair_R{wc.read}.fastq.gz")
    output:
        bam=temp(pjoin(OUTGROUP_MAP_TMP, "{sample_id}.unpaired_R{read}.bam"))
    threads: MAPPING_OUTGROUP_THREADS
    resources:
        mem_mb=OUTGROUP_MAP_UNPAIRED_MEM_MB,
        runtime=OUTGROUP_MAP_UNPAIRED_RUNTIME
    conda:
        "../envs/mapper.yaml"
    message:
        "Map unpaired R{wildcards.read} for outgroup {wildcards.sample_id} using {MAPPER_CMD}"
    shell:
        r"""
        mkdir -p {OUTGROUP_MAP_TMP}
        {MAPPER_CMD} -t {threads} {MAP_EXTRA} {MAPPER_REF_ARG} {input.fq} \
            | samtools view -Sb - \
            > {output.bam}
        """

def _merge_inputs(wc):
    bams = [pjoin(OUTGROUP_MAP_TMP, f"{wc.sample_id}.paired.bam")]
    if MAP_UNPAIRED:
        bams += [
            pjoin(OUTGROUP_MAP_TMP, f"{wc.sample_id}.unpaired_R1.bam"),
            pjoin(OUTGROUP_MAP_TMP, f"{wc.sample_id}.unpaired_R2.bam"),
        ]
    return bams

rule outgroup_final_bam:
    """
    Merge, sort, and index final outgroup BAMs from short-read alignments.
    """
    wildcard_constraints:
        sample_id=_wc_regex(SHORTREAD_SAMPLES)
    input:
        _merge_inputs
    output:
        bam=pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam"),
        bai=pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam.bai")
    threads: MAPPING_OUTGROUP_THREADS
    resources:
        mem_mb=OUTGROUP_MAP_FINAL_MEM_MB,
        runtime=OUTGROUP_MAP_FINAL_RUNTIME
    conda:
        "../envs/samtools.yaml"
    message:
        "Merge, sort, and index BAM for outgroup {wildcards.sample_id}"
    shell:
        r"""
        mkdir -p {OUTGROUP_FINAL_DIR}
        tmp=$(mktemp -p {OUTGROUP_MAP_TMP} {wildcards.sample_id}.merged.XXXXXX.bam)
        set -- {input}
        if [ "$#" -gt 1 ]; then
            samtools merge -f "$tmp" "$@"
        else
            cp "$1" "$tmp"
        fi
        samtools sort -@ {threads} -o {output.bam} "$tmp"
        samtools index {output.bam}
        rm -f "$tmp"
        """

#########################
# Long-read mapping flow #
#########################

def _lr_input_fastq(wc):
    # Prefer filtered fastq if filtering is enabled; otherwise use merged input
    if LR_CFG.get("filter", {}).get("enabled", True):
        return pjoin(OUTGROUP_LR_FILTER_DIR, f"{wc.sample_id}.fastq.gz")
    else:
        return f"{OUTGROUP_MERGED_DIR}/{wc.sample_id}.fastq.gz"

rule filtlong_outgroup:
    """
    Filter long-read outgroup FASTQs with filtlong before mapping.
    """
    wildcard_constraints:
        sample_id=_wc_regex(LONGREAD_SAMPLES)
    input:
        fq=lambda wc: f"{OUTGROUP_MERGED_DIR}/{wc.sample_id}.fastq.gz"
    output:
        fq=pjoin(OUTGROUP_LR_FILTER_DIR, "{sample_id}.fastq.gz")
    params:
        min_length=lambda wc: LR_CFG.get("filter", {}).get("min_length", 1000),
        keep_percent=lambda wc: LR_CFG.get("filter", {}).get("keep_percent", 90),
        extra=lambda wc: LR_CFG.get("filter", {}).get("extra", "")
    threads: OUTGROUP_FILTLONG_THREADS
    conda:
        "../envs/mapper.yaml"
    message:
        "Filter long reads (filtlong) for outgroup {wildcards.sample_id}"
    shell:
        r"""
        mkdir -p {OUTGROUP_LR_FILTER_DIR}
        filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {params.extra} \
            {input.fq} | gzip -c > {output.fq}
        """

rule map_outgroup_long_minimap2:
    """
    Map outgroup long reads with minimap2 using per-sample presets.
    """
    wildcard_constraints:
        sample_id=_wc_regex(LONGREAD_SAMPLES)
    input:
        ref=lambda wc: REF_UNZIPPED,
        idx=lambda wc: REF_MAP_INDEX_OUTGROUP,
        fq=_lr_input_fastq
    output:
        bam=pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam"),
        bai=pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam.bai")
    params:
        preset=lambda wc: OUTGROUP_MINIMAP2_PRESET.get(wc.sample_id, "map-pb"),
        extra=lambda wc: LR_CFG.get("mapping", {}).get("extra", "")
    threads: MAPPING_OUTGROUP_THREADS
    resources:
        mem_mb=OUTGROUP_LONGREAD_MAP_MEM_MB,
        runtime=OUTGROUP_LONGREAD_MAP_RUNTIME
    conda:
        "../envs/mapper.yaml"
    message:
        "Map long reads for outgroup {wildcards.sample_id} using minimap2 ({params.preset})"
    shell:
        r"""
        mkdir -p {OUTGROUP_FINAL_DIR}
        minimap2 -t {threads} -x {params.preset} {params.extra} {input.ref} {input.fq} \
            | samtools view -Sb - \
            | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

############################
# Aggregate outgroup BAMs  #
############################

rule outgroup_bams:
    """
    Aggregate: produce final BAMs for all short-read outgroup samples.
    """
    input:
        expand(pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam"), sample_id=[sid for sid in OUTGROUP_READ_TYPE.keys() if OUTGROUP_READ_TYPE[sid] in ("short","long")])

rule manifest_outgroup_longread_filter_storage:
    """
    Document the managed storage location for filtered outgroup long reads.
    """
    input:
        expand(pjoin(OUTGROUP_LR_FILTER_DIR, "{sample_id}.fastq.gz"), sample_id=LONGREAD_SAMPLES) if LR_CFG.get("filter", {}).get("enabled", True) else []
    output:
        readme=manifest_paths(OUTGROUP_LR_FILTER_DIR)[0],
        yaml=manifest_paths(OUTGROUP_LR_FILTER_DIR)[1]
    run:
        write_storage_manifest(
            OUTGROUP_LR_FILTER_DIR,
            "Outgroup Long-read Filter Storage",
            "manifest_outgroup_longread_filter_storage",
            {
                "asset_type": "filtered long-read outgroup FASTQ files",
                "sample_count": len(LONGREAD_SAMPLES),
                "filter_enabled": LR_CFG.get("filter", {}).get("enabled", True),
            },
        )

rule manifest_outgroup_bam_storage:
    """
    Document the managed storage location for final outgroup BAM files.
    """
    input:
        expand(pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam"), sample_id=[sid for sid in OUTGROUP_READ_TYPE.keys() if OUTGROUP_READ_TYPE[sid] in ("short", "long")]),
        expand(pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam.bai"), sample_id=[sid for sid in OUTGROUP_READ_TYPE.keys() if OUTGROUP_READ_TYPE[sid] in ("short", "long")])
    output:
        readme=manifest_paths(OUTGROUP_FINAL_DIR)[0],
        yaml=manifest_paths(OUTGROUP_FINAL_DIR)[1]
    run:
        write_storage_manifest(
            OUTGROUP_FINAL_DIR,
            "Outgroup BAM Storage",
            "manifest_outgroup_bam_storage",
            {
                "asset_type": "final outgroup BAM and BAI files",
                "sample_count": len(OUTGROUP_SAMPLE_IDS),
                "trim_dir": OUTGROUP_TRIM_DIR,
                "longread_filter_dir": OUTGROUP_LR_FILTER_DIR,
                "mapping_tmp_dir": OUTGROUP_MAP_TMP,
            },
        )


OUTGROUP_ALL_MAPPED_SAMPLES = [sid for sid in OUTGROUP_READ_TYPE.keys() if OUTGROUP_READ_TYPE[sid] in ("short", "long")]
OUTGROUP_MAPPING_TARGETS = [
    *expand(pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam"), sample_id=OUTGROUP_ALL_MAPPED_SAMPLES),
    *expand(pjoin(OUTGROUP_FINAL_DIR, "{sample_id}.bam.bai"), sample_id=OUTGROUP_ALL_MAPPED_SAMPLES),
    *rules.manifest_outgroup_bam_storage.output,
    *rules.manifest_outgroup_longread_filter_storage.output,
]
