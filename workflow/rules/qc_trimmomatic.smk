"""
QC, trimming (Trimmomatic), and post-trim QC for outgroup SRA reads.

Inputs:
- Paired FASTQs per outgroup sample from merge_fastq in outgroup.smk

Config keys used (with defaults shown):
- qc.outgroup.fastqc.contaminants: null
- qc.outgroup.fastqc.adapters: null
- qc.outgroup.trimmomatic.adapters_fa: null  # e.g., data/adapters/TruSeq3-PE.fa
- qc.outgroup.trimmomatic.clip: "2:30:10"
- qc.outgroup.trimmomatic.slidingwindow: "4:20"
- qc.outgroup.trimmomatic.leading: 3
- qc.outgroup.trimmomatic.trailing: 3
- qc.outgroup.trimmomatic.minlen: 36
- qc.outgroup.trimmomatic.extra: ""   # extra args appended to trimmomatic
"""

from os.path import join as pjoin
from pathlib import Path

# Resolve QC config with optional scope-specific overrides
def _qc_scope(scope: str):
    base = {k: v for k, v in (config.get("qc", {}) or {}).items() if k not in {"ingroup", "outgroup"}}
    scoped = (config.get("qc", {}) or {}).get(scope, {})
    if isinstance(scoped, dict):
        merged = base.copy()
        merged.update(scoped)
        return merged
    return base


def _fastqc_output_stem(path):
    name = Path(str(path)).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _outgroup_pre_fastqc_html(sample_id, read):
    return pjoin(PRE_QC_DIR, f"{_fastqc_output_stem(f'{OUTGROUP_MERGED_DIR}/{sample_id}_{read}.fastq.gz')}_fastqc.html")


def _outgroup_pre_fastqc_zip(sample_id, read):
    return pjoin(PRE_QC_DIR, f"{_fastqc_output_stem(f'{OUTGROUP_MERGED_DIR}/{sample_id}_{read}.fastq.gz')}_fastqc.zip")


def _outgroup_post_fastqc_html(sample_id, read):
    return pjoin(POST_QC_DIR, f"{sample_id}_pair_R{read}_fastqc.html")


def _outgroup_post_fastqc_zip(sample_id, read):
    return pjoin(POST_QC_DIR, f"{sample_id}_pair_R{read}_fastqc.zip")

OUTGROUP_QC_CFG = _qc_scope("outgroup")

# Only operate on outgroups detected as short-read in the shared sample helpers.
SHORTREAD_SAMPLES = [sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "short"]

PRE_QC_DIR = pjoin(OUTGROUP_QC_BASE_DIR, "pre")
POST_QC_DIR = pjoin(OUTGROUP_QC_BASE_DIR, "post")
OUTGROUP_TRIM_DIR_RULE = OUTGROUP_TRIM_DIR

FASTQC_CONTAM = (OUTGROUP_QC_CFG.get("fastqc", {}) or {}).get("contaminants", None)
FASTQC_ADAPTERS = (OUTGROUP_QC_CFG.get("fastqc", {}) or {}).get("adapters", None)
OUTGROUP_FASTQC_THREADS = _resolve_named_threads(
    4,
    "qc",
    "outgroup_fastqc",
    legacy_section=(OUTGROUP_QC_CFG.get("fastqc", {}) or {}),
    legacy_fallback_global=False,
)
TRIM_CFG = OUTGROUP_QC_CFG.get("trimmomatic", {}) or {}
OUTGROUP_TRIM_MEM_MB = _resolve_mem_mb(200000, "qc", "outgroup_trimmomatic")
OUTGROUP_TRIM_RUNTIME = _resolve_runtime(1440, "qc", "outgroup_trimmomatic")

TRIM_ADAPTERS = TRIM_CFG.get("adapters_fa", None)
TRIM_CLIP = TRIM_CFG.get("clip", "2:30:10")
TRIM_SLIDING = TRIM_CFG.get("slidingwindow", "4:20")
TRIM_LEADING = TRIM_CFG.get("leading", 3)
TRIM_TRAILING = TRIM_CFG.get("trailing", 3)
TRIM_MINLEN = TRIM_CFG.get("minlen", 36)
TRIM_EXTRA = str(TRIM_CFG.get("extra", "") or "").strip()

rule fastqc_outgroup_pre:
    """
    Run FastQC on merged outgroup reads before trimming.
    """
    input:
        fq=lambda wc: f"{OUTGROUP_MERGED_DIR}/{wc.sample_id}_{wc.read}.fastq.gz",
        contaminants=lambda wc: FASTQC_CONTAM if FASTQC_CONTAM else [],
        adapters=lambda wc: FASTQC_ADAPTERS if FASTQC_ADAPTERS else []
    output:
        html=pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.html"),
        zip=pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.zip")
    params:
        outdir=PRE_QC_DIR,
        contaminants=FASTQC_CONTAM,
        adapters=FASTQC_ADAPTERS,
        contam_opt=lambda wc: (f"--contaminants {FASTQC_CONTAM}" if FASTQC_CONTAM else ""),
        adapter_opt=lambda wc: (f"--adapters {FASTQC_ADAPTERS}" if FASTQC_ADAPTERS else ""),
        actual_html=lambda wc: _outgroup_pre_fastqc_html(wc.sample_id, wc.read),
        actual_zip=lambda wc: _outgroup_pre_fastqc_zip(wc.sample_id, wc.read)
    threads: OUTGROUP_FASTQC_THREADS
    conda:
        "../envs/fastqc.yaml"
    message:
        "FastQC (pre-trim) on outgroup {wildcards.sample_id} R{wildcards.read}"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} \
            {params.adapter_opt} \
            {params.contam_opt} \
            {input.fq}
        mv -f {params.actual_html} {output.html}
        mv -f {params.actual_zip} {output.zip}
        """

rule multiqc_outgroup_pre:
    """
    Aggregate pre-trim FastQC reports for short-read outgroups with MultiQC.
    """
    input:
        expand(pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.html"), sample_id=SHORTREAD_SAMPLES, read=["1","2"])
    output:
        html=pjoin(PRE_QC_DIR, "multiqc_report.html")
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        mkdir -p {PRE_QC_DIR}
        multiqc {PRE_QC_DIR} --outdir {PRE_QC_DIR} --filename multiqc_report.html
        """

rule trimmomatic_outgroup_pe:
    """
    Trim paired-end outgroup reads with Trimmomatic.
    """
    input:
        fq1=lambda wc: f"{OUTGROUP_MERGED_DIR}/{wc.sample_id}_1.fastq.gz",
        fq2=lambda wc: f"{OUTGROUP_MERGED_DIR}/{wc.sample_id}_2.fastq.gz",
        adapters=lambda wc: TRIM_ADAPTERS if TRIM_ADAPTERS else []
    output:
        pair1=pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"),
        unpair1=pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_unpair_R1.fastq.gz"),
        pair2=pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"),
        unpair2=pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_unpair_R2.fastq.gz")
    params:
        adapters=TRIM_ADAPTERS,
        clip=TRIM_CLIP,
        sliding=TRIM_SLIDING,
        leading=TRIM_LEADING,
        trailing=TRIM_TRAILING,
        minlen=TRIM_MINLEN,
        extra=TRIM_EXTRA,
        outdir=OUTGROUP_TRIM_DIR_RULE
    threads: QC_OUTGROUP_THREADS
    resources:
        mem_mb=OUTGROUP_TRIM_MEM_MB,
        runtime=OUTGROUP_TRIM_RUNTIME
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
    """
    Run FastQC on paired trimmed outgroup reads.
    """
    input:
        fq=lambda wc: pjoin(OUTGROUP_TRIM_DIR_RULE, f"{wc.sample_id}_pair_R{wc.read}.fastq.gz"),
        contaminants=lambda wc: FASTQC_CONTAM if FASTQC_CONTAM else [],
        adapters=lambda wc: FASTQC_ADAPTERS if FASTQC_ADAPTERS else []
    output:
        html=pjoin(POST_QC_DIR, "{sample_id}_R{read}_fastqc.html"),
        zip=pjoin(POST_QC_DIR, "{sample_id}_R{read}_fastqc.zip")
    params:
        outdir=POST_QC_DIR,
        contaminants=FASTQC_CONTAM,
        adapters=FASTQC_ADAPTERS,
        contam_opt=lambda wc: (f"--contaminants {FASTQC_CONTAM}" if FASTQC_CONTAM else ""),
        adapter_opt=lambda wc: (f"--adapters {FASTQC_ADAPTERS}" if FASTQC_ADAPTERS else ""),
        actual_html=lambda wc: _outgroup_post_fastqc_html(wc.sample_id, wc.read),
        actual_zip=lambda wc: _outgroup_post_fastqc_zip(wc.sample_id, wc.read)
    threads: OUTGROUP_FASTQC_THREADS
    conda:
        "../envs/fastqc.yaml"
    message:
        "FastQC (post-trim) on outgroup {wildcards.sample_id} R{wildcards.read}"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} \
            {params.adapter_opt} \
            {params.contam_opt} \
            {input.fq}
        mv -f {params.actual_html} {output.html}
        mv -f {params.actual_zip} {output.zip}
        """

rule multiqc_outgroup_post:
    """
    Aggregate post-trim FastQC reports for short-read outgroups with MultiQC.
    """
    input:
        expand(pjoin(POST_QC_DIR, "{sample_id}_R{read}_fastqc.html"), sample_id=SHORTREAD_SAMPLES, read=["1","2"])
    output:
        html=pjoin(POST_QC_DIR, "multiqc_report.html")
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        mkdir -p {POST_QC_DIR}
        multiqc {POST_QC_DIR} --outdir {POST_QC_DIR} --filename multiqc_report.html
        """

rule outgroup_qc_trim_all:
    """
    Aggregate: run pre-QC, trimming, and post-QC for all short-read outgroup samples.
    """
    input:
        # pre-QC htmls for R1/R2
        expand(pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.html"), sample_id=SHORTREAD_SAMPLES, read=["1","2"]),
        # trimmed pairs
        expand(pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"), sample_id=SHORTREAD_SAMPLES),
        expand(pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"), sample_id=SHORTREAD_SAMPLES),
        # post-QC multiqc
        pjoin(POST_QC_DIR, "multiqc_report.html")

rule manifest_outgroup_trim_storage:
    """
    Document the managed storage location for trimmed outgroup FASTQs.
    """
    input:
        expand(pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"), sample_id=SHORTREAD_SAMPLES),
        expand(pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"), sample_id=SHORTREAD_SAMPLES)
    output:
        readme=manifest_paths(OUTGROUP_TRIM_DIR_RULE)[0],
        yaml=manifest_paths(OUTGROUP_TRIM_DIR_RULE)[1]
    run:
        write_storage_manifest(
            OUTGROUP_TRIM_DIR_RULE,
            "Outgroup Trimmed Read Storage",
            "manifest_outgroup_trim_storage",
            {
                "asset_type": "trimmed outgroup FASTQ files",
                "sample_count": len(SHORTREAD_SAMPLES),
                "source_dir": OUTGROUP_MERGED_DIR,
            },
        )

rule manifest_outgroup_qc_storage:
    """
    Document the managed storage location for outgroup QC reports.
    """
    input:
        pre_multiqc=pjoin(PRE_QC_DIR, "multiqc_report.html"),
        post_multiqc=pjoin(POST_QC_DIR, "multiqc_report.html")
    output:
        readme=manifest_paths(OUTGROUP_QC_BASE_DIR)[0],
        yaml=manifest_paths(OUTGROUP_QC_BASE_DIR)[1]
    run:
        write_storage_manifest(
            OUTGROUP_QC_BASE_DIR,
            "Outgroup QC Storage",
            "manifest_outgroup_qc_storage",
            {
                "asset_type": "FastQC and MultiQC outputs for outgroup reads",
                "pre_qc_dir": PRE_QC_DIR,
                "post_qc_dir": POST_QC_DIR,
                "sample_count": len(SHORTREAD_SAMPLES),
            },
        )


OUTGROUP_QC_TARGETS = [
    *expand(pjoin(PRE_QC_DIR, "{sample_id}_{read}_fastqc.html"), sample_id=SHORTREAD_SAMPLES, read=["1", "2"]),
    *expand(pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"), sample_id=SHORTREAD_SAMPLES),
    *expand(pjoin(OUTGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"), sample_id=SHORTREAD_SAMPLES),
    pjoin(PRE_QC_DIR, "multiqc_report.html"),
    pjoin(POST_QC_DIR, "multiqc_report.html"),
    *rules.manifest_outgroup_trim_storage.output,
    *rules.manifest_outgroup_qc_storage.output,
] if SHORTREAD_SAMPLES else []
