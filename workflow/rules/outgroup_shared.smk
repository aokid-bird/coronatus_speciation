"""
Shared outgroup-derived metadata rules.

These rules consume sample-level outgroup FASTQs regardless of whether they were
created locally or provided externally, so they can be loaded in both local and
cluster environments.

Config keys used:
- outgroups
- storage.outgroup.merged_dir
"""

LONGREAD_SAMPLES = [sid for sid, typ in OUTGROUP_READ_TYPE.items() if typ == "long"]


rule manifest_outgroup_merged_storage:
    """
    Record provenance for merged sample-level outgroup FASTQ files.
    """
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
