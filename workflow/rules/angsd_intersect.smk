# rules/angsd.smk

def _bams_for_group(group):
    return [
        f"{config_bam_dir}/{s}.bam"
        for s in _SAMPLES_DF[_SAMPLES_DF[group_col] == group]["sample"].astype(str)
    ]


# Rules related to angsd_intersect
rule make_bamlist_group:
    """
    Generate bamlists for each groups defined
    """
    input:
        samples = config["samples"],
        bams = lambda wc: _bams_for_group(wc.group)
    output:
        bamlist = f"results/bamlists/{output_prefix}/{{group}}/bamlist.txt"
    run:
        df = _SAMPLES_DF
        selected = df[df[group_col] == wildcards.group]
        bams = selected["sample"].apply(lambda s: f"{config_bam_dir}/{s}.bam")
        bams.to_csv(output.bamlist, index=False, header=False)

rule angsd_intersect_group:
    """
    Population-specific ANGSD to find intersecting sites among populations 
    """
    input:
        bamlist = f"results/bamlists/{output_prefix}/{{group}}/bamlist.txt"
    output:
        geno = f"results/angsd_intersect/{output_prefix}/{{group}}/gl.geno.gz"
    log:
        f"logs/{output_prefix}/angsd_intersect_{{group}}.log"
    params:
        outprefix = lambda wc: f"results/angsd_intersect/{output_prefix}/{wc.group}/gl",
        refpath = REF,
        minInd = lambda wc: max(1, int(
            sum(1 for _ in open(f"results/bamlists/{output_prefix}/{wc.group}/bamlist.txt")) * config['minIndRatio']['intersect']
        )),
        extra=config["angsd_common_args"].strip() + " " + config["angsd_args"]["intersect"].strip()
    threads: config["threads"]
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        angsd -out {params.outprefix} -b {input.bamlist} \
              -ref {params.refpath} -anc {params.refpath} \
              -minInd {params.minInd} \
              {params.extra} \
              -nThreads {threads} 2> {log}
        """

rule intersect_sites_group:
    input:
        genos = expand(f"results/angsd_intersect/{output_prefix}/{{group}}/gl.geno.gz", group=groups)
    output:
        sites = f"results/intersect_sites/{output_prefix}/intersect.txt",
        scafs = f"results/intersect_sites/{output_prefix}/intersect.chr"
    params:
        script="workflow/scripts/intersect_sites.py"
    shell:
        """
        conda run -n snakemake_intersect_sites python {params.script} \
            --genofiles {input.genos} \
            --out_sites {output.sites} \
            --out_scafs {output.scafs}
        """

rule angsd_sites_index:
    input:
        sites = f"results/intersect_sites/{output_prefix}/intersect.txt"
    output:
        sites_idx = f"results/intersect_sites/{output_prefix}/intersect.txt.bin"
    log:
        f"logs/{output_prefix}/angsd_intersect_sites_index.log"
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        angsd sites index {input.sites}
        """
