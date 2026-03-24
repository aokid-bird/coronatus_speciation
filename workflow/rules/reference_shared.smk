# rules/reference_shared.smk

"""
Shared reference-derived outputs and mapper indices.

These rules are environment-agnostic and may be used in either local or cluster
runs once the reference FASTA path has been resolved.
"""

rule derive_reference_chroms:
    """
    Generate a chromosome or contig list from the FASTA index for ANGSD site filters.
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
    """
    Export contig lengths from the FASTA index into a CSV summary for downstream analyses.
    """
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


_reference_manifest_inputs = [
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


rule manifest_reference_storage:
    """
    Record provenance for the resolved reference storage location and derived index files.
    """
    input:
        _reference_manifest_inputs
    output:
        readme=manifest_paths(REFERENCE_MANIFEST_DIR)[0],
        yaml=manifest_paths(REFERENCE_MANIFEST_DIR)[1]
    run:
        write_storage_manifest(
            REFERENCE_MANIFEST_DIR,
            "Reference Storage",
            "manifest_reference_storage",
            {
                "asset_type": "reference FASTA and mapper index files",
                "reference_directory": REFERENCE_DIR,
                "reference_fasta": REF,
                "reference_fai": FAI,
                "reference_chr": REF_CHR,
                "reference_source": REFERENCE_ACCESSION or REF_CONFIG_PATH or "user-supplied",
                "download_method": REFERENCE_METHOD,
            },
        )


from pathlib import Path as _Path
if _Path(REF).suffix == ".gz":
    rule index_reference_samtools:
        """
        Create samtools faidx indices for a gzipped FASTA reference.
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

    rule decompress_reference_fasta:
        """
        Create an uncompressed FASTA copy when mappers require one.
        """
        input:
            fasta_gz=REF
        output:
            fasta=REF_UNZIPPED
        shell:
            r"""
            gzip -dc {input.fasta_gz} > {output.fasta}
            """


rule index_bwa:
    """
    Build bwa index files for the resolved reference FASTA.
    """
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
        runtime=480
    shell:
        r"""
        mkdir -p logs/reference $(dirname {output})
        export OMP_NUM_THREADS={threads}
        bwa index -p {REF_INDEX_PREFIX_BWA} {input.fasta} &> {log}
        """


rule index_bwa_mem2:
    """
    Build bwa-mem2 index files for the resolved reference FASTA.
    """
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
        runtime=480
    shell:
        r"""
        mkdir -p logs/reference $(dirname {output})
        export OMP_NUM_THREADS={threads}
        bwa-mem2 index -p {REF_INDEX_PREFIX_BWAMEM2} {input.fasta} &> {log}
        """
