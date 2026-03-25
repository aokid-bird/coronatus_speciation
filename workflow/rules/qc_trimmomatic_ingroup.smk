"""
QC, trimming (Trimmomatic), and post-trim QC for ingroup reads listed in config.samples.

Inputs:
- Paired FASTQs per ingroup sample resolved from config.samples metadata
  using per-sample directory/prefix/suffix/extension fields with legacy fallback
  to reads.ingroup_dir and {sample}_1/{sample}_2.fastq.gz naming.

Config keys (scoped under qc.ingroup with fallback to qc base):
- qc.ingroup.fastqc.contaminants
- qc.ingroup.fastqc.adapters
- qc.ingroup.trimmomatic.adapters_fa
- qc.ingroup.trimmomatic.clip
- qc.ingroup.trimmomatic.slidingwindow
- qc.ingroup.trimmomatic.leading
- qc.ingroup.trimmomatic.trailing
- qc.ingroup.trimmomatic.minlen
- qc.ingroup.trimmomatic.extra
"""

from os.path import join as pjoin
from pathlib import Path
import re


def _qc_scope(scope: str):
    base = {k: v for k, v in (config.get("qc", {}) or {}).items() if k not in {"ingroup", "outgroup"}}
    scoped = (config.get("qc", {}) or {}).get(scope, {})
    if isinstance(scoped, dict):
        merged = base.copy()
        merged.update(scoped)
        return merged
    return base


def _wc_regex(ids):
    return "(" + "|".join(map(re.escape, ids)) + ")" if ids else r"a^"


def _fastqc_output_stem(path):
    name = Path(str(path)).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _ingroup_pre_fastqc_html(sample_id, read):
    return pjoin(PRE_QC_DIR, f"{_fastqc_output_stem(ingroup_fastq_path(sample_id, read))}_fastqc.html")


def _ingroup_pre_fastqc_zip(sample_id, read):
    return pjoin(PRE_QC_DIR, f"{_fastqc_output_stem(ingroup_fastq_path(sample_id, read))}_fastqc.zip")


def _ingroup_post_fastqc_html(sample_id, read):
    return pjoin(POST_QC_DIR, f"{sample_id}_pair_R{read}_fastqc.html")


def _ingroup_post_fastqc_zip(sample_id, read):
    return pjoin(POST_QC_DIR, f"{sample_id}_pair_R{read}_fastqc.zip")


INGROUP_QC_CFG = _qc_scope("ingroup")

PRE_QC_DIR = pjoin(INGROUP_QC_BASE_DIR, "pre")
POST_QC_DIR = pjoin(INGROUP_QC_BASE_DIR, "post")
INGROUP_TRIM_DIR_RULE = INGROUP_TRIM_DIR

FASTQC_CONTAM = (INGROUP_QC_CFG.get("fastqc", {}) or {}).get("contaminants", None)
FASTQC_ADAPTERS = (INGROUP_QC_CFG.get("fastqc", {}) or {}).get("adapters", None)
INGROUP_FASTQC_THREADS = _resolve_named_threads(
    4,
    "qc",
    "ingroup_fastqc",
    legacy_section=(INGROUP_QC_CFG.get("fastqc", {}) or {}),
    legacy_fallback_global=False,
)
TRIM_CFG = INGROUP_QC_CFG.get("trimmomatic", {}) or {}
INGROUP_TRIM_MEM_MB = _resolve_mem_mb(200000, "qc", "ingroup_trimmomatic")
INGROUP_TRIM_RUNTIME = _resolve_runtime(1440, "qc", "ingroup_trimmomatic")

TRIM_ADAPTERS = TRIM_CFG.get("adapters_fa", None)
TRIM_CLIP = TRIM_CFG.get("clip", "2:30:10")
TRIM_SLIDING = TRIM_CFG.get("slidingwindow", "4:20")
TRIM_LEADING = TRIM_CFG.get("leading", 3)
TRIM_TRAILING = TRIM_CFG.get("trailing", 3)
TRIM_MINLEN = TRIM_CFG.get("minlen", 36)
TRIM_EXTRA = str(TRIM_CFG.get("extra", "") or "").strip()


rule fastqc_ingroup_pre:
    """
    Run FastQC on ingroup reads before trimming.
    """
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_STORAGE_SAMPLE_IDS),
        read="1|2"
    input:
        fq=lambda wc: ingroup_fastq_path(wc.sample_id, wc.read),
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
        actual_html=lambda wc: _ingroup_pre_fastqc_html(wc.sample_id, wc.read),
        actual_zip=lambda wc: _ingroup_pre_fastqc_zip(wc.sample_id, wc.read)
    threads: INGROUP_FASTQC_THREADS
    conda:
        "../envs/fastqc.yaml"
    message:
        "FastQC (pre-trim) on ingroup {wildcards.sample_id} R{wildcards.read}"
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


rule multiqc_ingroup_pre:
    """
    Aggregate pre-trim ingroup FastQC reports with MultiQC.
    """
    input:
        expand(_ingroup_pre_fastqc_html("{sample_id}", "{read}"), sample_id=INGROUP_STORAGE_SAMPLE_IDS, read=["1","2"])
    output:
        html=pjoin(PRE_QC_DIR, "multiqc_report.html")
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        multiqc {PRE_QC_DIR} --outdir {PRE_QC_DIR}
        """


rule trimmomatic_ingroup_pe:
    """
    Trim paired-end ingroup reads with Trimmomatic.
    """
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_STORAGE_SAMPLE_IDS)
    input:
        fq1=lambda wc: ingroup_fastq_path(wc.sample_id, "1"),
        fq2=lambda wc: ingroup_fastq_path(wc.sample_id, "2"),
        adapters=lambda wc: TRIM_ADAPTERS if TRIM_ADAPTERS else []
    output:
        pair1=pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"),
        unpair1=pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_unpair_R1.fastq.gz"),
        pair2=pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"),
        unpair2=pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_unpair_R2.fastq.gz")
    params:
        adapters=TRIM_ADAPTERS,
        clip=TRIM_CLIP,
        sliding=TRIM_SLIDING,
        leading=TRIM_LEADING,
        trailing=TRIM_TRAILING,
        minlen=TRIM_MINLEN,
        extra=TRIM_EXTRA,
        outdir=INGROUP_TRIM_DIR_RULE
    threads: QC_INGROUP_THREADS
    resources:
        mem_mb=INGROUP_TRIM_MEM_MB,
        runtime=INGROUP_TRIM_RUNTIME
    conda:
        "../envs/trimmomatic.yaml"
    message:
        "Trimmomatic PE on ingroup {wildcards.sample_id}"
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


rule fastqc_ingroup_post:
    """
    Run FastQC on paired trimmed ingroup reads.
    """
    wildcard_constraints:
        sample_id=_wc_regex(INGROUP_STORAGE_SAMPLE_IDS),
        read="1|2"
    input:
        fq=lambda wc: pjoin(INGROUP_TRIM_DIR_RULE, f"{wc.sample_id}_pair_R{wc.read}.fastq.gz"),
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
        actual_html=lambda wc: _ingroup_post_fastqc_html(wc.sample_id, wc.read),
        actual_zip=lambda wc: _ingroup_post_fastqc_zip(wc.sample_id, wc.read)
    threads: INGROUP_FASTQC_THREADS
    conda:
        "../envs/fastqc.yaml"
    message:
        "FastQC (post-trim) on ingroup {wildcards.sample_id} R{wildcards.read}"
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


rule multiqc_ingroup_post:
    """
    Aggregate post-trim ingroup FastQC reports with MultiQC.
    """
    input:
        expand(_ingroup_post_fastqc_html("{sample_id}", "{read}"), sample_id=INGROUP_STORAGE_SAMPLE_IDS, read=["1","2"])
    output:
        html=pjoin(POST_QC_DIR, "multiqc_report.html")
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        multiqc {POST_QC_DIR} --outdir {POST_QC_DIR}
        """


rule ingroup_qc_trim_all:
    """
    Aggregate: run pre-QC, trimming, and post-QC for all ingroup samples.
    """
    input:
        # pre-QC htmls for R1/R2
        expand(_ingroup_pre_fastqc_html("{sample_id}", "{read}"), sample_id=INGROUP_STORAGE_SAMPLE_IDS, read=["1","2"]),
        # trimmed pairs
        expand(pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"), sample_id=INGROUP_STORAGE_SAMPLE_IDS),
        expand(pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"), sample_id=INGROUP_STORAGE_SAMPLE_IDS),
        # post-QC multiqc
        pjoin(POST_QC_DIR, "multiqc_report.html")

rule manifest_ingroup_trim_storage:
    """
    Document the managed storage location for trimmed ingroup FASTQs.
    """
    input:
        expand(pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"), sample_id=INGROUP_STORAGE_SAMPLE_IDS),
        expand(pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"), sample_id=INGROUP_STORAGE_SAMPLE_IDS)
    output:
        readme=manifest_paths(INGROUP_TRIM_DIR_RULE)[0],
        yaml=manifest_paths(INGROUP_TRIM_DIR_RULE)[1]
    run:
        write_storage_manifest(
            INGROUP_TRIM_DIR_RULE,
            "Ingroup Trimmed Read Storage",
            "manifest_ingroup_trim_storage",
            {
                "asset_type": "trimmed ingroup FASTQ files",
                "sample_count": len(INGROUP_STORAGE_SAMPLE_IDS),
                "source_metadata_tsv": config["samples"],
                "legacy_default_dir": INGROUP_READS_DIR,
            },
        )

rule manifest_ingroup_qc_storage:
    """
    Document the managed storage location for ingroup QC reports.
    """
    input:
        pre_multiqc=pjoin(PRE_QC_DIR, "multiqc_report.html"),
        post_multiqc=pjoin(POST_QC_DIR, "multiqc_report.html")
    output:
        readme=manifest_paths(INGROUP_QC_BASE_DIR)[0],
        yaml=manifest_paths(INGROUP_QC_BASE_DIR)[1]
    run:
        write_storage_manifest(
            INGROUP_QC_BASE_DIR,
            "Ingroup QC Storage",
            "manifest_ingroup_qc_storage",
            {
                "asset_type": "FastQC and MultiQC outputs for ingroup reads",
                "pre_qc_dir": PRE_QC_DIR,
                "post_qc_dir": POST_QC_DIR,
                "sample_count": len(INGROUP_STORAGE_SAMPLE_IDS),
            },
        )


INGROUP_QC_TARGETS = [
    *expand(_ingroup_pre_fastqc_html("{sample_id}", "{read}"), sample_id=INGROUP_STORAGE_SAMPLE_IDS, read=["1", "2"]),
    *expand(pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R1.fastq.gz"), sample_id=INGROUP_STORAGE_SAMPLE_IDS),
    *expand(pjoin(INGROUP_TRIM_DIR_RULE, "{sample_id}_pair_R2.fastq.gz"), sample_id=INGROUP_STORAGE_SAMPLE_IDS),
    pjoin(PRE_QC_DIR, "multiqc_report.html"),
    pjoin(POST_QC_DIR, "multiqc_report.html"),
    *rules.manifest_ingroup_trim_storage.output,
    *rules.manifest_ingroup_qc_storage.output,
]
