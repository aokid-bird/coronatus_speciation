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


rule samples_with_outgroups_metadata:
    """
    Extend the sample metadata with selected outgroups for downstream plots.
    """
    input:
        samples=config["samples"],
        outgroups=config["outgroups"]
    output:
        samples_aug=f"results/metadata/{output_prefix}/samples_plus_outgroups.tsv"
    run:
        import os
        import pandas as pd

        os.makedirs(os.path.dirname(output.samples_aug), exist_ok=True)
        samples_df = pd.read_csv(input.samples, sep="\t")
        outgroups_df = pd.read_csv(input.outgroups, sep="\t")

        cols = samples_df.columns.tolist()
        if "sample" not in cols:
            raise ValueError("samples TSV must include a 'sample' column")

        rows = []
        for _, row in outgroups_df.iterrows():
            new_row = {col: None for col in cols}
            new_row["sample"] = str(row["sample_id"]) if "sample_id" in outgroups_df.columns else None
            if group_col in cols and "taxon" in outgroups_df.columns:
                new_row[group_col] = str(row["taxon"]) if pd.notna(row["taxon"]) else None
            rows.append(new_row)

        augmented = pd.concat([samples_df, pd.DataFrame(rows)], ignore_index=True)
        augmented.to_csv(output.samples_aug, sep="\t", index=False)
