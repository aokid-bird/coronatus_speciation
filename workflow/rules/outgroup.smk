"""
rules/outgroup.smk
Fetch outgroup SRA reads, handling paired short-reads and single-end long-reads separately.
"""

import re as _re

def _wc_rgx(ids):
    return "(" + "|".join(map(_re.escape, ids)) + ")" if ids else r"a^"

# Build SRR lists by read type from OUTGROUP_SAMPLES and OUTGROUP_READ_TYPE (from common.smk)
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
    Fetch paired-end short-read FASTQs from SRA and gzip
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
        if [ "{config[environment]}" != "local" ]; then
            echo "This rule must be run in local environment." >&2
            exit 1
        fi

        if [ -s {output.fq1} ] && [ -s {output.fq2} ]; then
            echo "[{wildcards.srr}] Already exists." > {log}
        else
            # Choose compressor (prefer pigz for speed)
            if command -v pigz >/dev/null 2>&1; then
                COMPRESSOR="pigz -p {threads} -f"
            else
                COMPRESSOR="gzip -f"
            fi

            # If uncompressed FASTQs already exist, skip re-download and just gzip
            mkdir -p {OUTGROUP_RAW_DIR}
            if [ -s "{OUTGROUP_RAW_DIR}/{wildcards.srr}_1.fastq" ] && [ -s "{OUTGROUP_RAW_DIR}/{wildcards.srr}_2.fastq" ]; then
                echo "[{wildcards.srr}] Found existing FASTQ; skipping download and gzipping." >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_1.fastq" >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_2.fastq" >> {log} 2>&1
            else
                # Download FASTQ (uncompressed) and split mates
                fasterq-dump {wildcards.srr} --split-files -e {threads} -O {OUTGROUP_RAW_DIR} &> {log}
                # Compress paired-end files
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_1.fastq" >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}_2.fastq" >> {log} 2>&1
            fi

            # Verify expected outputs for paired-end datasets
            if [ ! -s {output.fq1} ] || [ ! -s {output.fq2} ]; then
                echo "[{wildcards.srr}] ERROR: Expected paired .fastq.gz not found after compression." >> {log} 2>&1
                exit 1
            fi
        fi
        """


rule fetch_sra_single:
    """
    Fetch single-end (long-read) FASTQ from SRA and gzip
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
        if [ "{config[environment]}" != "local" ]; then
            echo "This rule must be run in local environment." >&2
            exit 1
        fi

        if [ -s {output.fq} ]; then
            echo "[{wildcards.srr}] Already exists." > {log}
        else
            # Choose compressor (prefer pigz for speed)
            if command -v pigz >/dev/null 2>&1; then
                COMPRESSOR="pigz -p {threads} -f"
            else
                COMPRESSOR="gzip -f"
            fi

            # If uncompressed FASTQ already exists, skip re-download and just gzip
            mkdir -p {OUTGROUP_RAW_DIR}
            if [ -s "{OUTGROUP_RAW_DIR}/{wildcards.srr}.fastq" ]; then
                echo "[{wildcards.srr}] Found existing FASTQ; skipping download and gzipping." >> {log} 2>&1
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}.fastq" >> {log} 2>&1
            else
                # Download FASTQ (uncompressed) without splitting
                fasterq-dump {wildcards.srr} -e {threads} -O {OUTGROUP_RAW_DIR} &> {log}
                $COMPRESSOR "{OUTGROUP_RAW_DIR}/{wildcards.srr}.fastq" >> {log} 2>&1
            fi

            # Verify expected output exists
            if [ ! -s {output.fq} ]; then
                echo "[{wildcards.srr}] ERROR: Expected single .fastq.gz not found after compression." >> {log} 2>&1
                exit 1
            fi
        fi
        """

rule merge_fastq:
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

# Optional: merge single-end long-read fastqs per sample (if present under {OUTGROUP_RAW_DIR}/{srr}.fastq.gz)
LONGREAD_SAMPLES = [sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "long"]

rule merge_fastq_long:
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

rule manifest_outgroup_merged_storage:
    input:
        shortread_merged_1=expand(f"{OUTGROUP_MERGED_DIR}/{{sample_id}}_1.fastq.gz", sample_id=[sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "short"]),
        shortread_merged_2=expand(f"{OUTGROUP_MERGED_DIR}/{{sample_id}}_2.fastq.gz", sample_id=[sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "short"]),
        longread_merged=expand(f"{OUTGROUP_MERGED_DIR}/{{sample_id}}.fastq.gz", sample_id=LONGREAD_SAMPLES)
    output:
        readme=manifest_paths(OUTGROUP_MERGED_DIR)[0],
        yaml=manifest_paths(OUTGROUP_MERGED_DIR)[1]
    run:
        write_storage_manifest(
            OUTGROUP_MERGED_DIR,
            "Outgroup Merged FASTQ Storage",
            "manifest_outgroup_merged_storage",
            {
                "asset_type": "merged outgroup FASTQ files",
                "sample_count": len(OUTGROUP_SAMPLE_IDS),
                "shortread_samples": len([sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "short"]),
                "longread_samples": len(LONGREAD_SAMPLES),
            },
        )
