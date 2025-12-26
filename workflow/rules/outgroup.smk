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
        fq1="data/raw/outgroup/{srr}_1.fastq.gz",
        fq2="data/raw/outgroup/{srr}_2.fastq.gz"
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
            if [ -s "data/raw/outgroup/{wildcards.srr}_1.fastq" ] && [ -s "data/raw/outgroup/{wildcards.srr}_2.fastq" ]; then
                echo "[{wildcards.srr}] Found existing FASTQ; skipping download and gzipping." >> {log} 2>&1
                $COMPRESSOR "data/raw/outgroup/{wildcards.srr}_1.fastq" >> {log} 2>&1
                $COMPRESSOR "data/raw/outgroup/{wildcards.srr}_2.fastq" >> {log} 2>&1
            else
                # Download FASTQ (uncompressed) and split mates
                fasterq-dump {wildcards.srr} --split-files -e {threads} -O data/raw/outgroup/ &> {log}
                # Compress paired-end files
                $COMPRESSOR "data/raw/outgroup/{wildcards.srr}_1.fastq" >> {log} 2>&1
                $COMPRESSOR "data/raw/outgroup/{wildcards.srr}_2.fastq" >> {log} 2>&1
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
        fq="data/raw/outgroup/{srr}.fastq.gz"
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
            if [ -s "data/raw/outgroup/{wildcards.srr}.fastq" ]; then
                echo "[{wildcards.srr}] Found existing FASTQ; skipping download and gzipping." >> {log} 2>&1
                $COMPRESSOR "data/raw/outgroup/{wildcards.srr}.fastq" >> {log} 2>&1
            else
                # Download FASTQ (uncompressed) without splitting
                fasterq-dump {wildcards.srr} -e {threads} -O data/raw/outgroup/ &> {log}
                $COMPRESSOR "data/raw/outgroup/{wildcards.srr}.fastq" >> {log} 2>&1
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
        fq1s=lambda wc: expand("data/raw/outgroup/{srr}_1.fastq.gz", srr=OUTGROUP_SAMPLES[wc.sample_id]),
        fq2s=lambda wc: expand("data/raw/outgroup/{srr}_2.fastq.gz", srr=OUTGROUP_SAMPLES[wc.sample_id])
    output:
        fq1="data/merged/outgroup/{sample_id}_1.fastq.gz",
        fq2="data/merged/outgroup/{sample_id}_2.fastq.gz"
    shell:
        """
        cat {input.fq1s} > {output.fq1}
        cat {input.fq2s} > {output.fq2}
        """

# Optional: merge single-end long-read fastqs per sample (if present under data/raw/outgroup/{srr}.fastq.gz)
LONGREAD_SAMPLES = [sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "long"]

rule merge_fastq_long:
    wildcard_constraints:
        sample_id=_wc_rgx(LONGREAD_SAMPLES)
    input:
        fqs=lambda wc: expand("data/raw/outgroup/{srr}.fastq.gz", srr=OUTGROUP_SAMPLES[wc.sample_id])
    output:
        fq="data/merged/outgroup/{sample_id}.fastq.gz"
    shell:
        """
        cat {input.fqs} > {output.fq}
        """
