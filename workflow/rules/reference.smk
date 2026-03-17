# rules/reference.smk

"""
Reference acquisition and indexing.

Local-only download/extract rules are gated by REF_FROM_METADATA, but generic
indexing rules (samtools faidx, bwa/bwa-mem2 index) are always available so
cluster runs can generate required indices when provided a reference path.

Outputs:
 - {REF} (gzipped FASTA) [local mode]
 - {FAI}, {GZI} (samtools faidx indices)
 - {REF_CHR} (contig list for ANGSD -rf convenience)
 - Mapper index for bwa or bwa-mem2 (selected via REF_MAP_INDEX_MAIN)
"""

if REF_FROM_METADATA:
    rule download_reference_zip:
        """
        Download reference archive via NCBI Datasets, or wget if URL provided.
        """
        output:
            zip_file=REFERENCE_ZIP
        log:
            f"logs/fetch_reference/{REFERENCE_ACCESSION}.log"
        conda:
            "../envs/ncbi_datasets.yaml"
        threads: 2
        message:
            f"Fetching reference {REFERENCE_ACCESSION} via {REFERENCE_METHOD}"
        shell:
            r"""
            mkdir -p $(dirname {output.zip_file}) logs/fetch_reference
            if [ -n "{REFERENCE_URL}" ]; then
                wget -O {output.zip_file} "{REFERENCE_URL}" &> {log}
            else
                datasets download genome accession {REFERENCE_ACCESSION} \
                  --include genome \
                  --filename {output.zip_file} \
                  &> {log}
            fi
            """

    rule extract_reference_fna:
        """
        Extract FASTA from the zip and concatenate to a single gzipped FASTA.
        """
        input:
            zip_file=REFERENCE_ZIP
        output:
            fasta_gz=REF
        conda:
            "../envs/ncbi_datasets.yaml"
        shell:
            r"""
            tmpdir=$(mktemp -d)
            unzip -q -o -d "$tmpdir" {input.zip_file}
            # Concatenate all .fna files in the archive
            find "$tmpdir" -type f -name "*.fna" -print0 | xargs -0 cat | gzip -c > {output.fasta_gz}
            rm -rf "$tmpdir"
            """

    rule index_reference_samtools:
        """
        Create samtools faidx indices (.fai and .gzi for gzipped FASTA).
        """
        input:
            fasta_gz=REF
        output:
            fai=FAI,
            gzi=GZI if GZI is not None else temp(".ref.gzi.notused")
        conda:
            "../envs/samtools.yaml"
        threads: 2
        shell:
            r"""
            tmp=$(mktemp -t refbgzip).bgz
            gunzip -c {input.fasta_gz} | bgzip -c > "$tmp"
            mv "$tmp" {input.fasta_gz}
            samtools faidx {input.fasta_gz}
            """

rule derive_reference_chroms:
    """
    Generate a chromosome/contig list from .fai (for ANGSD -rf convenience).
    """
    input:
        fai=FAI
    output:
        chr=REF_CHR
    shell:
        """
        cut -f1 {input.fai} > {output.chr}
        """


rule export_contig_lengths:
    """Export contig lengths derived from the FASTA index (.fai)."""
    input:
        fai=FAI
    output:
        csv=f"results/reference/{output_prefix}/contig_lengths.csv"
    params:
        script="workflow/scripts/export_contig_lengths.py"
    conda:
        "../envs/python_pandas.yaml"
    shell:
        """
        python {params.script} --fai {input.fai} --output {output.csv}
        """

_prepare_reference_inputs = [
    path
    for path in (
        REFERENCE_ZIP if REF_FROM_METADATA and REFERENCE_ZIP else None,
        REF,
        FAI,
        REF_CHR,
        REF_UNZIPPED,
        *REF_MAP_INDEXES,
    )
    if path
]


rule prepare_reference:
    """
    Aggregate rule ensuring reference is downloaded and indexed (local use).
    """
    input:
        _prepare_reference_inputs
    output:
        # fasta=REF,
        # fai=FAI,
        # chr=REF_CHR,
        # mapindex=REF_MAP_INDEX_MAIN
    message:
        "Reference prepared (downloaded + samtools indexed + contigs listed + mapper index)"

# Generic indexing rules (always available)

# Ensure we have an uncompressed FASTA for mapping indices if the provided
# reference is gzipped. Define this rule only when REF is gz to avoid
# input==output edge cases when REF is already uncompressed.
from pathlib import Path as _Path
if _Path(REF).suffix == ".gz":
    rule decompress_reference_fasta:
        input:
            fasta_gz=REF
        output:
            fasta=REF_UNZIPPED
        shell:
            r"""
            gzip -dc {input.fasta_gz} > {output.fasta}
            """

rule index_bwa:
    input:
        fasta=REF_UNZIPPED
    output:
        BWA_INDEX_MAIN
    log:
        "logs/reference/index_bwa.log"
    conda:
        "../envs/mapper.yaml"
    threads: 2
    resources:
        mem_mb=64000,
        runtime="08:00:00"
    shell:
        r"""
        mkdir -p logs/reference $(dirname {output})
        export OMP_NUM_THREADS={threads}
        bwa index -p {REF_INDEX_PREFIX_BWA} {input.fasta} &> {log}
        """

rule index_bwa_mem2:
    input:
        fasta=REF_UNZIPPED
    output:
        BWAMEM2_INDEX_MAIN
    log:
        "logs/reference/index_bwa_mem2.log"
    conda:
        "../envs/mapper.yaml"
    threads: 2
    resources:
        mem_mb=64000,
        runtime="08:00:00"
    shell:
        r"""
        mkdir -p logs/reference $(dirname {output})
        export OMP_NUM_THREADS={threads}
        # Index to a writable prefix if configured
        bwa-mem2 index -p {REF_INDEX_PREFIX_BWAMEM2} {input.fasta} &> {log}
        """
