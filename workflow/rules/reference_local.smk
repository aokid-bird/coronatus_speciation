"""
Local-only reference acquisition.

These rules handle downloading and extracting the reference FASTA when the
config does not provide a ready-to-use path. Shared indexing and derived
reference outputs are defined in reference_shared.smk so cluster runs do not
load local-only download behavior.

Config keys used:
- reference.fasta
- references_tsv
- storage.reference.dir
- reference_download_method
"""

if REF_FROM_METADATA:
    rule download_reference_zip:
        """
        Download the reference archive locally via NCBI Datasets or wget.
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
        Extract FASTA records from the downloaded archive into the project reference path.
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
            find "$tmpdir" -type f -name "*.fna" -print0 | xargs -0 cat | gzip -c > {output.fasta_gz}
            rm -rf "$tmpdir"
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
    Ensure the local reference has been acquired and all shared derived reference files exist.
    """
    input:
        _prepare_reference_inputs
    message:
        "Reference prepared locally"
