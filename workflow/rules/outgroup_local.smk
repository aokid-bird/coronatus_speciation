"""
Local-only outgroup FASTQ acquisition and merge rules.

This module is intentionally limited to local bootstrap work: fetching SRA
FASTQs, merging SRR-level files into sample-level files, and writing storage
manifests for those local assets.
"""

import re as _re


def _wc_rgx(ids):
    return "(" + "|".join(map(_re.escape, ids)) + ")" if ids else r"a^"


SHORTREAD_SRR_IDS = [
    srr
    for sid, typ in OUTGROUP_READ_TYPE.items()
    if typ == "short"
    for srr in OUTGROUP_SAMPLES[sid]
]
LONGREAD_SRR_IDS = [
    srr
    for sid, typ in OUTGROUP_READ_TYPE.items()
    if typ == "long"
    for srr in OUTGROUP_SAMPLES[sid]
]


rule fetch_sra_paired:
    """
    Download paired-end outgroup reads from SRA into the local raw outgroup directory.
    """
    wildcard_constraints:
        srr=_wc_rgx(SHORTREAD_SRR_IDS)
    output:
        fq1=f"{OUTGROUP_RAW_DIR}/{{srr}}_1.fastq.gz",
        fq2=f"{OUTGROUP_RAW_DIR}/{{srr}}_2.fastq.gz"
    log:
        "logs/fetch_sra/{srr}.log"
    conda:
        "../envs/sra_tools.yaml"
    threads: 4
    message:
        "Downloading SRA paired-end {wildcards.srr}"
    shell:
        """
        if [ -s {output.fq1} ] && [ -s {output.fq2} ]; then
            echo "[{wildcards.srr}] Already exists." > {log}
        else
            if command -v pigz >/dev/null 2>&1; then
                COMPRESSOR="pigz -p {threads} -f"
            else
                COMPRESSOR="gzip -f"
            fi

            mkdir -p {OUTGROUP_RAW_DIR}
            if [ -s "{OUTGROUP_RAW_DIR}/{wildcards.srr}_1.fastq" ] && [ -s "{OUTGROUP_RAW_DIR}/{wildcards.srr}_2.fastq" ]; then
                echo "[{wildcards.srr}] Found existing FASTQ; skipping download and gzipping." >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_1.fastq" >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_2.fastq" >> {log} 2>&1
            else
                fasterq-dump {wildcards.srr} --split-files -e {threads} -O {OUTGROUP_RAW_DIR} &> {log}
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_1.fastq" >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_2.fastq" >> {log} 2>&1
            fi

            if [ ! -s {output.fq1} ] || [ ! -s {output.fq2} ]; then
                echo "[{wildcards.srr}] ERROR: Expected paired .fastq.gz not found after compression." >> {log} 2>&1
                exit 1
            fi
        fi
        """


rule fetch_sra_single:
    """
    Download single-end outgroup reads from SRA into the local raw outgroup directory.
    """
    wildcard_constraints:
        srr=_wc_rgx(LONGREAD_SRR_IDS)
    output:
        fq=f"{OUTGROUP_RAW_DIR}/{{srr}}.fastq.gz"
    log:
        "logs/fetch_sra/{srr}.log"
    conda:
        "../envs/sra_tools.yaml"
    threads: 4
    message:
        "Downloading SRA single-end {wildcards.srr}"
    shell:
        """
        if [ -s {output.fq} ]; then
            echo "[{wildcards.srr}] Already exists." > {log}
        else
            if command -v pigz >/dev/null 2>&1; then
                COMPRESSOR="pigz -p {threads} -f"
            else
                COMPRESSOR="gzip -f"
            fi

            mkdir -p {OUTGROUP_RAW_DIR}
            if [ -s "{OUTGROUP_RAW_DIR}/{wildcards.srr}.fastq" ]; then
                echo "[{wildcards.srr}] Found existing FASTQ; skipping download and gzipping." >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}.fastq" >> {log} 2>&1
            else
                fasterq-dump {wildcards.srr} -e {threads} -O {OUTGROUP_RAW_DIR} &> {log}
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}.fastq" >> {log} 2>&1
            fi

            if [ ! -s {output.fq} ]; then
                echo "[{wildcards.srr}] ERROR: Expected single .fastq.gz not found after compression." >> {log} 2>&1
                exit 1
            fi
        fi
        """


rule merge_fastq:
    """
    Merge paired-end outgroup FASTQs across SRR runs into one sample-level FASTQ pair.
    """
    input:
        fq1s=lambda wc: [f"{OUTGROUP_RAW_DIR}/{srr}_1.fastq.gz" for srr in OUTGROUP_SAMPLES[wc.sample_id]],
        fq2s=lambda wc: [f"{OUTGROUP_RAW_DIR}/{srr}_2.fastq.gz" for srr in OUTGROUP_SAMPLES[wc.sample_id]]
    output:
        fq1=f"{OUTGROUP_MERGED_DIR}/{{sample_id}}_1.fastq.gz",
        fq2=f"{OUTGROUP_MERGED_DIR}/{{sample_id}}_2.fastq.gz"
    shell:
        """
        mkdir -p {OUTGROUP_MERGED_DIR}
        cat {input.fq1s} > {output.fq1}
        cat {input.fq2s} > {output.fq2}
        """


rule merge_fastq_long:
    """
    Merge single-end long-read outgroup FASTQs across SRR runs into one sample-level FASTQ.
    """
    wildcard_constraints:
        sample_id=_wc_rgx(LONGREAD_SAMPLES)
    input:
        fqs=lambda wc: [f"{OUTGROUP_RAW_DIR}/{srr}.fastq.gz" for srr in OUTGROUP_SAMPLES[wc.sample_id]]
    output:
        fq=f"{OUTGROUP_MERGED_DIR}/{{sample_id}}.fastq.gz"
    shell:
        """
        mkdir -p {OUTGROUP_MERGED_DIR}
        cat {input.fqs} > {output.fq}
        """


rule manifest_outgroup_raw_storage:
    """
    Record provenance for locally downloaded outgroup raw FASTQ files.
    """
    input:
        paired=expand(f"{OUTGROUP_RAW_DIR}/{{srr}}_1.fastq.gz", srr=SHORTREAD_SRR_IDS),
        paired_mates=expand(f"{OUTGROUP_RAW_DIR}/{{srr}}_2.fastq.gz", srr=SHORTREAD_SRR_IDS),
        single=expand(f"{OUTGROUP_RAW_DIR}/{{srr}}.fastq.gz", srr=LONGREAD_SRR_IDS)
    output:
        readme=manifest_paths(OUTGROUP_RAW_DIR)[0],
        yaml=manifest_paths(OUTGROUP_RAW_DIR)[1]
    run:
        write_storage_manifest(
            OUTGROUP_RAW_DIR,
            "Outgroup Raw FASTQ Storage",
            "manifest_outgroup_raw_storage",
            {
                "asset_type": "downloaded outgroup SRA FASTQ files",
                "sample_count": len(OUTGROUP_SAMPLE_IDS),
                "srr_count": len(OUTGROUP_SRR_IDS),
                "read_types": ",".join(sorted(set(OUTGROUP_READ_TYPE.values()))),
            },
        )
