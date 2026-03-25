"""
Mapping trimmed ingroup reads to the reference.

Inputs:
- Trimmed ingroup FASTQs from qc_trimmomatic_ingroup.smk
- Shared reference FASTA and mapper indices from reference_shared.smk

Config keys used:
- mapping.mapper
- mapping.ingroup.mapper
- mapping.ingroup.map_unpaired
- mapping.ingroup.extra
- storage.shared.bam.ingroup_dir
- storage.bam.ingroup_dir
- resources.mapping.ingroup_paired
- resources.mapping.ingroup_unpaired
- resources.mapping.ingroup_merge
- threads.mapping.ingroup

Outputs final coordinate-sorted BAM and index in the resolved ingroup BAM storage directory.
"""

import re
from os.path import join as pjoin

TRIM_DIR = INGROUP_TRIM_DIR
MAP_TMP = INGROUP_MAP_TMP_DIR
FINAL_DIR = INGROUP_BAM_STORAGE_DIR
ALIAS_DIR = config_bam_dir

MAP_CFG = MAPPING_INGROUP_CFG
MAP_UNPAIRED = bool(MAP_CFG.get("map_unpaired", True))
MAP_EXTRA = str(MAP_CFG.get("extra", "") or "").strip()
INGROUP_MAP_PAIRED_MEM_MB = _resolve_mem_mb(64000, "mapping", "ingroup_paired")
INGROUP_MAP_PAIRED_RUNTIME = _resolve_runtime(1440, "mapping", "ingroup_paired")
INGROUP_MAP_UNPAIRED_MEM_MB = _resolve_mem_mb(64000, "mapping", "ingroup_unpaired")
INGROUP_MAP_UNPAIRED_RUNTIME = _resolve_runtime(1440, "mapping", "ingroup_unpaired")
INGROUP_MAP_FINAL_MEM_MB = _resolve_mem_mb(200000, "mapping", "ingroup_merge")
INGROUP_MAP_FINAL_RUNTIME = _resolve_runtime(1440, "mapping", "ingroup_merge")

MAPPER_CMD = "bwa-mem2 mem" if MAPPER_INGROUP == "bwa-mem2" else "bwa mem"
MAPPER_REF_ARG = REF_MAP_ARG_INGROUP


def _wc_regex(ids):
    return "(" + "|".join(map(re.escape, ids)) + ")" if ids else r"a^"


INGROUP_ANALYSIS_BAM_TARGETS = (
    expand(pjoin(ALIAS_DIR, "{sample_id}.bam"), sample_id=INGROUP_SAMPLE_IDS)
    if INGROUP_BAM_ALIAS_ACTIVE
    else expand(pjoin(FINAL_DIR, "{sample_id}.bam"), sample_id=INGROUP_STORAGE_SAMPLE_IDS)
)
INGROUP_ANALYSIS_BAI_TARGETS = (
    expand(pjoin(ALIAS_DIR, "{sample_id}.bam.bai"), sample_id=INGROUP_SAMPLE_IDS)
    if INGROUP_BAM_ALIAS_ACTIVE
    else expand(pjoin(FINAL_DIR, "{sample_id}.bam.bai"), sample_id=INGROUP_STORAGE_SAMPLE_IDS)
)


rule map_ingroup_paired:
    """
    Map paired trimmed ingroup reads with the configured short-read mapper.
    """
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_STORAGE_SAMPLE_IDS)
    input:
        ref=lambda wc: REF_UNZIPPED,
        idx=lambda wc: REF_MAP_INDEX_INGROUP,
        fq1=lambda wc: pjoin(TRIM_DIR, f"{wc.sample_id}_pair_R1.fastq.gz"),
        fq2=lambda wc: pjoin(TRIM_DIR, f"{wc.sample_id}_pair_R2.fastq.gz")
    output:
        bam=temp(pjoin(MAP_TMP, "{sample_id}.paired.bam"))
    threads: MAPPING_INGROUP_THREADS
    resources:
        mem_mb=INGROUP_MAP_PAIRED_MEM_MB,
        runtime=INGROUP_MAP_PAIRED_RUNTIME
    conda:
        "../envs/mapper.yaml"
    message:
        "Map paired reads for ingroup {wildcards.sample_id} using {MAPPER_CMD}"
    shell:
        r"""
        mkdir -p {MAP_TMP}
        {MAPPER_CMD} -t {threads} {MAP_EXTRA} {MAPPER_REF_ARG} {input.fq1} {input.fq2} \
            | samtools view -Sb - \
            > {output.bam}
        """


rule map_ingroup_unpaired:
    """
    Map trimmed ingroup unpaired reads when unpaired mapping is enabled.
    """
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_STORAGE_SAMPLE_IDS)
    input:
        ref=lambda wc: REF_UNZIPPED,
        idx=lambda wc: REF_MAP_INDEX_INGROUP,
        fq=lambda wc: pjoin(TRIM_DIR, f"{wc.sample_id}_unpair_R{wc.read}.fastq.gz")
    output:
        bam=temp(pjoin(MAP_TMP, "{sample_id}.unpaired_R{read}.bam"))
    threads: MAPPING_INGROUP_THREADS
    resources:
        mem_mb=INGROUP_MAP_UNPAIRED_MEM_MB,
        runtime=INGROUP_MAP_UNPAIRED_RUNTIME
    conda:
        "../envs/mapper.yaml"
    message:
        "Map unpaired R{wildcards.read} for ingroup {wildcards.sample_id} using {MAPPER_CMD}"
    shell:
        r"""
        mkdir -p {MAP_TMP}
        {MAPPER_CMD} -t {threads} {MAP_EXTRA} {MAPPER_REF_ARG} {input.fq} \
            | samtools view -Sb - \
            > {output.bam}
        """


def _merge_inputs(wc):
    bams = [pjoin(MAP_TMP, f"{wc.sample_id}.paired.bam")]
    if MAP_UNPAIRED:
        bams += [
            pjoin(MAP_TMP, f"{wc.sample_id}.unpaired_R1.bam"),
            pjoin(MAP_TMP, f"{wc.sample_id}.unpaired_R2.bam"),
        ]
    return bams


rule ingroup_final_bam:
    """
    Merge, sort, and index final ingroup BAMs from short-read alignments.
    """
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_STORAGE_SAMPLE_IDS)
    input:
        _merge_inputs
    output:
        bam=pjoin(FINAL_DIR, "{sample_id}.bam"),
        bai=pjoin(FINAL_DIR, "{sample_id}.bam.bai")
    threads: MAPPING_INGROUP_THREADS
    resources:
        mem_mb=INGROUP_MAP_FINAL_MEM_MB,
        runtime=INGROUP_MAP_FINAL_RUNTIME
    conda:
        "../envs/samtools.yaml"
    message:
        "Merge, sort, and index BAM for ingroup {wildcards.sample_id}"
    shell:
        r"""
        mkdir -p {FINAL_DIR}
        tmp=$(mktemp -p {MAP_TMP} {wildcards.sample_id}.merged.XXXXXX.bam)
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


rule ingroup_bams:
    """
    Aggregate: produce final BAMs for all ingroup samples.
    """
    input:
        INGROUP_ANALYSIS_BAM_TARGETS


if INGROUP_BAM_ALIAS_ACTIVE:
    rule ingroup_bam_alias:
        """
        Expose user-facing ingroup BAM aliases for downstream analyses.
        """
        wildcard_constraints:
            sample_id=_wc_regex(INGROUP_SAMPLE_IDS)
        input:
            bam=lambda wc: ingroup_storage_bam_path(ingroup_storage_sample_id(wc.sample_id)),
            bai=lambda wc: f"{ingroup_storage_bam_path(ingroup_storage_sample_id(wc.sample_id))}.bai"
        output:
            bam=pjoin(ALIAS_DIR, "{sample_id}.bam"),
            bai=pjoin(ALIAS_DIR, "{sample_id}.bam.bai")
        message:
            "Link ingroup BAM alias for downstream sample {wildcards.sample_id}"
        shell:
            r"""
            mkdir -p {ALIAS_DIR}
            ln -sfn {input.bam} {output.bam}
            ln -sfn {input.bai} {output.bai}
            """

rule manifest_ingroup_bam_storage:
    """
    Document the managed storage location for final ingroup BAM files.
    """
    input:
        expand(pjoin(FINAL_DIR, "{sample_id}.bam"), sample_id=INGROUP_STORAGE_SAMPLE_IDS),
        expand(pjoin(FINAL_DIR, "{sample_id}.bam.bai"), sample_id=INGROUP_STORAGE_SAMPLE_IDS)
    output:
        readme=manifest_paths(FINAL_DIR)[0],
        yaml=manifest_paths(FINAL_DIR)[1]
    run:
        write_storage_manifest(
            FINAL_DIR,
            "Ingroup BAM Storage",
            "manifest_ingroup_bam_storage",
            {
                "asset_type": "final ingroup BAM and BAI files",
                "sample_count": len(INGROUP_STORAGE_SAMPLE_IDS),
                "trim_dir": TRIM_DIR,
                "mapping_tmp_dir": MAP_TMP,
            },
        )


INGROUP_MAPPING_TARGETS = [
    *INGROUP_ANALYSIS_BAM_TARGETS,
    *INGROUP_ANALYSIS_BAI_TARGETS,
    *rules.manifest_ingroup_bam_storage.output,
]
