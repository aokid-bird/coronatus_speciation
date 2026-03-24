"""
Build union sites across ingroup populations from ANGSD intersect outputs,
then slice outgroup BAMs to those sites to reduce size for downstream ANGSD.
"""

rule union_sites_groups:
    input:
        genos = expand(f"results/angsd_intersect/{output_prefix}/{{group}}/gl.geno.gz", group=groups)
    output:
        sites = f"results/union_sites/{output_prefix}/union.sites",
        scafs = f"results/union_sites/{output_prefix}/union.chr",
        bed   = f"results/union_sites/{output_prefix}/union.bed"
    params:
        script = "workflow/scripts/union_sites.py"
    conda:
        "../envs/intersect_sites.yaml"
    shell:
        r"""
        python {params.script} --genofiles {input.genos} --out_sites {output.sites} --out_scafs {output.scafs}
        # TSV (1-based pos) -> BED (0-based start, 1-based end) single-base intervals
        awk 'BEGIN{{OFS="\t"}} {{print $1, $2-1, $2}}' {output.sites} > {output.bed}
        """


rule slice_outgroup_bam:
    """
    Slice each outgroup BAM to the union sites bed and index the result.
    """
    input:
        bam = f"{OUTGROUP_BAM_DIR}/{{sample_id}}.bam",
        bed = rules.union_sites_groups.output.bed
    output:
        bam = f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam",
        bai = f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam.bai"
    threads: SLICE_OUTGROUPS_THREADS
    resources:
        mem_mb=200000,
        runtime=1440
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {OUTGROUP_SLICED_DIR}
        # view reads overlapping union sites, sort, and index
        samtools view -@ {threads} -L {input.bed} -b {input.bam} \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule outgroup_bams_sliced:
    """
    Aggregate: ensure sliced BAMs are produced for all outgroup samples.
    """
    input:
        expand(f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam", sample_id=ACTIVE_OUTGROUP_SAMPLE_IDS)
