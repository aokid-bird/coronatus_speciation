"""
Mapping trimmed ingroup reads to the reference.

Outputs final coordinate-sorted BAM and index in config['bam_dir']/{sample}.bam.
"""

import re
from os.path import join as pjoin

TRIM_DIR = "results/trimmomatic/ingroup"
MAP_TMP = f"results/mapping/{output_prefix}/ingroup"
FINAL_DIR = config_bam_dir

MAP_CFG = MAPPING_INGROUP_CFG
MAP_UNPAIRED = bool(MAP_CFG.get("map_unpaired", True))
MAP_EXTRA = str(MAP_CFG.get("extra", "") or "").strip()

MAPPER_CMD = "bwa-mem2 mem" if MAPPER_INGROUP == "bwa-mem2" else "bwa mem"
MAPPER_REF_ARG = REF_MAP_ARG_INGROUP


def _wc_regex(ids):
    return "(" + "|".join(map(re.escape, ids)) + ")" if ids else r"a^"


rule map_ingroup_paired:
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_SAMPLE_IDS)
    input:
        ref=lambda wc: REF_UNZIPPED,
        idx=lambda wc: REF_MAP_INDEX_INGROUP,
        fq1=lambda wc: pjoin(TRIM_DIR, f"{wc.sample_id}_pair_R1.fastq.gz"),
        fq2=lambda wc: pjoin(TRIM_DIR, f"{wc.sample_id}_pair_R2.fastq.gz")
    output:
        bam=temp(pjoin(MAP_TMP, "{sample_id}.paired.bam"))
    threads: config.get("threads", 6)
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
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_SAMPLE_IDS)
    input:
        ref=lambda wc: REF_UNZIPPED,
        idx=lambda wc: REF_MAP_INDEX_INGROUP,
        fq=lambda wc: pjoin(TRIM_DIR, f"{wc.sample_id}_unpair_R{wc.read}.fastq.gz")
    output:
        bam=temp(pjoin(MAP_TMP, "{sample_id}.unpaired_R{read}.bam"))
    threads: config.get("threads", 6)
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
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_SAMPLE_IDS)
    input:
        _merge_inputs
    output:
        bam=pjoin(FINAL_DIR, "{sample_id}.bam"),
        bai=pjoin(FINAL_DIR, "{sample_id}.bam.bai")
    threads: config.get("threads", 6)
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
        expand(pjoin(FINAL_DIR, "{sample_id}.bam"), sample_id=INGROUP_SAMPLE_IDS)
