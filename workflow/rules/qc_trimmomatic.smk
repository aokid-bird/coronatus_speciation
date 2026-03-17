"""
QC, trimming (Trimmomatic), and post-trim QC for outgroup SRA reads.

Inputs:
- Paired FASTQs per outgroup sample from merge_fastq in outgroup.smk

Config keys used (with defaults shown):
- qc.outgroup.fastqc.contaminants: null
- qc.outgroup.trimmomatic.adapters_fa: null  # e.g., data/adapters/TruSeq3-PE.fa
- qc.outgroup.trimmomatic.clip: "2:30:10"
- qc.outgroup.trimmomatic.slidingwindow: "4:20"
- qc.outgroup.trimmomatic.leading: 3
- qc.outgroup.trimmomatic.trailing: 3
- qc.outgroup.trimmomatic.minlen: 36
- qc.outgroup.trimmomatic.extra: ""   # extra args appended to trimmomatic
"""

from os.path import join as pjoin

# Resolve QC config with optional scope-specific overrides
def _qc_scope(scope: str):
    base = {k: v for k, v in (config.get("qc", {}) or {}).items() if k not in {"ingroup", "outgroup"}}
    scoped = (config.get("qc", {}) or {}).get(scope, {})
    if isinstance(scoped, dict):
        merged = base.copy()
        merged.update(scoped)
        return merged
    return base

OUTGROUP_QC_CFG = _qc_scope("outgroup")

# Only operate on outgroups detected as short-read in common.smk
SHORTREAD_SAMPLES = [sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "short"]

PRE_QC_DIR = "results/qc/outgroup/pre"
POST_QC_DIR = "results/qc/outgroup/post"
TRIM_DIR = "results/trimmomatic/outgroup"

FASTQC_CONTAM = (OUTGROUP_QC_CFG.get("fastqc", {}) or {}).get("contaminants", None)
TRIM_CFG = OUTGROUP_QC_CFG.get("trimmomatic", {}) or {}

TRIM_ADAPTERS = TRIM_CFG.get("adapters_fa", None)
TRIM_CLIP = TRIM_CFG.get("clip", "2:30:10")
TRIM_SLIDING = TRIM_CFG.get("slidingwindow", "4:20")
TRIM_LEADING = TRIM_CFG.get("leading", 3)
TRIM_TRAILING = TRIM_CFG.get("trailing", 3)
TRIM_MINLEN = TRIM_CFG.get("minlen", 36)
TRIM_EXTRA = str(TRIM_CFG.get("extra", "") or "").strip()

rule fastqc_outgroup_pre:
    input:
        fq=lambda wc: f"data/merged/outgroup/{wc.sample_id}_{wc.read}.fastq.gz"
    output:
        html=pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.html"),
        zip=pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.zip")
    params:
        outdir=PRE_QC_DIR,
        contaminants=FASTQC_CONTAM,
        contam_opt=lambda wc: (f"--contaminants {FASTQC_CONTAM}" if FASTQC_CONTAM else "")
    threads: 4
    conda:
        "../envs/fastqc.yaml"
    message:
        "FastQC (pre-trim) on outgroup {wildcards.sample_id} R{wildcards.read}"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} \
            {params.contam_opt} \
            {input.fq}
        """

rule multiqc_outgroup_pre:
    input:
        expand(pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.html"), sample_id=SHORTREAD_SAMPLES, read=["1","2"])
    output:
        html=pjoin(PRE_QC_DIR, "multiqc_report.html")
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        multiqc {PRE_QC_DIR} --outdir {PRE_QC_DIR}
        """

rule trimmomatic_outgroup_pe:
    input:
        fq1=lambda wc: f"data/merged/outgroup/{wc.sample_id}_1.fastq.gz",
        fq2=lambda wc: f"data/merged/outgroup/{wc.sample_id}_2.fastq.gz"
    output:
        pair1=pjoin(TRIM_DIR, "{sample_id}_pair_R1.fastq.gz"),
        unpair1=pjoin(TRIM_DIR, "{sample_id}_unpair_R1.fastq.gz"),
        pair2=pjoin(TRIM_DIR, "{sample_id}_pair_R2.fastq.gz"),
        unpair2=pjoin(TRIM_DIR, "{sample_id}_unpair_R2.fastq.gz")
    params:
        adapters=TRIM_ADAPTERS,
        clip=TRIM_CLIP,
        sliding=TRIM_SLIDING,
        leading=TRIM_LEADING,
        trailing=TRIM_TRAILING,
        minlen=TRIM_MINLEN,
        extra=TRIM_EXTRA,
        outdir=TRIM_DIR
    threads: config.get("threads", 4)
    resources:
        mem_mb=200000,
        runtime="24:00:00"
    conda:
        "../envs/trimmomatic.yaml"
    message:
        "Trimmomatic PE on outgroup {wildcards.sample_id}"
    shell:
        r"""
        mkdir -p {params.outdir}
        if [ -n "{params.adapters}" ]; then
            ILLUM="ILLUMINACLIP:{params.adapters}:{params.clip}"
        else
            ILLUM=""
        fi
        trimmomatic PE -threads {threads} \
            {input.fq1} {input.fq2} \
            {output.pair1} {output.unpair1} \
            {output.pair2} {output.unpair2} \
            $ILLUM SLIDINGWINDOW:{params.sliding} \
            LEADING:{params.leading} TRAILING:{params.trailing} MINLEN:{params.minlen} {params.extra}
        """

rule fastqc_outgroup_post:
    input:
        fq=lambda wc: pjoin(TRIM_DIR, f"{wc.sample_id}_pair_R{wc.read}.fastq.gz")
    output:
        html=pjoin(POST_QC_DIR, "{sample_id}_R{read}_fastqc.html"),
        zip=pjoin(POST_QC_DIR, "{sample_id}_R{read}_fastqc.zip")
    params:
        outdir=POST_QC_DIR,
        contaminants=FASTQC_CONTAM,
        contam_opt=lambda wc: (f"--contaminants {FASTQC_CONTAM}" if FASTQC_CONTAM else "")
    threads: 4
    conda:
        "../envs/fastqc.yaml"
    message:
        "FastQC (post-trim) on outgroup {wildcards.sample_id} R{wildcards.read}"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} \
            {params.contam_opt} \
            {input.fq}
        """

rule multiqc_outgroup_post:
    input:
        expand(pjoin(POST_QC_DIR, "{sample_id}_R{read}_fastqc.html"), sample_id=SHORTREAD_SAMPLES, read=["1","2"])
    output:
        html=pjoin(POST_QC_DIR, "multiqc_report.html")
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        multiqc {POST_QC_DIR} --outdir {POST_QC_DIR}
        """

rule outgroup_qc_trim_all:
    """
    Aggregate: run pre-QC, trimming, and post-QC for all short-read outgroup samples.
    """
    input:
        # pre-QC htmls for R1/R2
        expand(pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.html"), sample_id=SHORTREAD_SAMPLES, read=["1","2"]),
        # trimmed pairs
        expand(pjoin(TRIM_DIR, "{sample_id}_pair_R1.fastq.gz"), sample_id=SHORTREAD_SAMPLES),
        expand(pjoin(TRIM_DIR, "{sample_id}_pair_R2.fastq.gz"), sample_id=SHORTREAD_SAMPLES),
        # post-QC multiqc
        pjoin(POST_QC_DIR, "multiqc_report.html")
