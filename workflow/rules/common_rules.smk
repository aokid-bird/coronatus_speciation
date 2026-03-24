"""
Small shared rules reused by more than one analysis block.
"""

rule bcf2vcf:
    """
    Convert an ANGSD BCF output into a compressed VCF for downstream tools.
    """
    input:
        bcf=lambda wc: f"results/angsd_{wc.angsd_run}/{output_prefix}/gl.bcf"
    output:
        vcf=f"results/angsd_{{angsd_run}}/{output_prefix}/gl.vcf.gz"
    conda:
        "../envs/bcftools_env.yaml"
    shell:
        """
        bcftools convert -O z -o {output.vcf} {input.bcf}
        """
