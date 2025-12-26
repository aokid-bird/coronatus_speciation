# scripts/make_bamlist.py
import pandas as pd

df = pd.read_csv(snakemake.input.samples, sep="\t")
selected = df[df[snakemake.params.group_col] == snakemake.wildcards.group]
bams = selected["sample"].apply(lambda s: f"{snakemake.params.bam_dir}/{s}.bam")
bams.to_csv(snakemake.output.bamlist, index=False, header=False)
