# Codex Integration Notes: QC, Trimming, and Mapping (Outgroups)

This document records the user prompts and Codex responses, plus a concise description of the changes made to integrate legacy R workflows into the Snakemake pipeline for outgroup QC, trimming, and mapping.

## User Prompts and Codex Responses

### Prompt 1 (User)
I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake.

0) General
Descriptions & Specifications:
- These codes were developed in different projects. Please adapt them to the present snakemake pipelines on config parameters, tsv files to be read, directories, etc.
- cluster.script in the old script intends to convert shell commands to PBS cluster. Please ignore this step since snakemake is intended to be ran in PBS cluster in most cases (except when environment == "local").
- They are ran under environment == "cluster"
- You may add columns to any already-present tsv files to retrieve metadata 
- Please adjust styles to snakemake. For example, loop inside R is not adequate but single rules with expand comand would be sufficient in this pipeline.

1) qc_trimmomatic.R
Descriptions:
This quality checks raeds, trims reads, and re-quality checks trimmed reads using FastQC, MultiQC, and Trimmomatic. 
Specifications:
- Ignore L. 17-54
- Parameters are set in the config file
- adapter files are specified in the config file

2) reference_mapping.R
Descriptions:
This script maps trimmed reads on reference fasta. This tries to map both paired (bwap) and unpaired (bwa1/bwa2).
Specifications:
- Use either bwa mem or bwa-mem2 depending on the config file specification
- Let users to specify parameters in the config file
- Let users to specify whether unpaired reads will be mapped and combined to paired mapped reads.

3) sra_manipulation.R
Descriptions:
This script intends to 0) download SRA, 1) trim SRA using Trimmomatic, 2) QCs trimmed reads, (1+2 for short reads), 3) trim SRA and QCs trimmed reads for PacBio long-reads, 4) map trimmed reads against reference, index them and combine them as single bams
Specifications:
- SRR fetch has been implemented, so please ignore
- trimming and bwa mapping specification could be switched depending on sequencer types defined in outgroups.tsv (the sequencer columns or any columns you think adequate)

### Response 1 (Codex)
I’ll scan the repo to understand the current Snakemake layout and configs, then draft targeted rules and config changes.

Summary of work completed in this step:
- Added rules for outgroup FastQC/MultiQC + Trimmomatic (short reads)
- Added mapping rules (bwa/bwa-mem2) including optional unpaired mapping and final BAM merge/sort/index
- Extended `common.smk` to infer outgroup read type (short vs long) based on `outgroups.tsv`
- Wired new rules into `Snakefile` and added minimal conda env files and config keys

### Prompt 2 (User)
Please do or think about the following
- bwa.yaml and mapper.yaml seem to be redundant. Please consolidate them and reflect changes in rules that use these environments
- add minimap2-based trimming/mapping for long reads and integrate those outputs as well

### Response 2 (Codex)
I’ll consolidate the mapping envs by switching reference indexing to use a single `mapper.yaml`, remove the old `bwa*.yaml` files, then add minimap2+filtlong and long-read rules (filter + minimap2 mapping) integrated into the outgroup BAM outputs.

Summary of work completed in this step:
- Consolidated mapping envs into `workflow/envs/mapper.yaml` (bwa, bwa-mem2, samtools, minimap2, filtlong)
- Removed `workflow/envs/bwa.yaml` and `workflow/envs/bwa_mem2.yaml`
- Switched reference indexing rules to use the consolidated env
- Added long-read flow: optional `filtlong` filtering and `minimap2` mapping using presets inferred from `sequencer`
- Aggregated final outgroup BAMs now include both short- and long-read samples

---

## Summary of Changes
- New short-read QC and trimming rules for outgroups using FastQC, MultiQC, and Trimmomatic.
- New short-read mapping rules using `bwa mem` or `bwa-mem2 mem`, with optional unpaired mapping and merged final BAMs.
- Consolidated mapping environment into a single conda env (`mapper.yaml`).
- Added long-read (PacBio/ONT) branch: optional read filtering with `filtlong` and mapping with `minimap2` using presets based on `sequencer` metadata.
- Extended `common.smk` to infer read type (short/long) and minimap2 presets from `data/outgroups.tsv`.
- Added/updated config keys for QC, Trimmomatic, mapping, and long-read options.

## Files Added or Updated
- New rules
  - `workflow/rules/qc_trimmomatic.smk`
    - `fastqc_outgroup_pre`, `multiqc_outgroup_pre`
    - `trimmomatic_outgroup_pe`
    - `fastqc_outgroup_post`, `multiqc_outgroup_post`
    - `outgroup_qc_trim_all` (aggregate)
  - `workflow/rules/mapping_outgroup.smk`
    - Short-read: `map_outgroup_paired`, `map_outgroup_unpaired`, `outgroup_final_bam`
    - Long-read: `filtlong_outgroup`, `map_outgroup_long_minimap2`
    - Aggregate: `outgroup_bams`
- Existing rules updated
  - `workflow/rules/common.smk`
    - Adds `OUTGROUP_READ_TYPE` (short vs long) from `outgroups.tsv` column `sequencer` (configurable)
    - Adds `OUTGROUP_MINIMAP2_PRESET` inference (`map-ont`, `map-hifi`, `map-pb`)
  - `workflow/rules/outgroup.smk`
    - Adds `merge_fastq_long` to concatenate single-end long-read fastqs
  - `workflow/rules/reference.smk`
    - Reference index rules switched to `../envs/mapper.yaml`
  - `workflow/Snakefile`
    - Includes new rule files in cluster mode
- Conda envs
  - Added: `workflow/envs/fastqc.yaml`, `workflow/envs/multiqc.yaml`, `workflow/envs/trimmomatic.yaml`, `workflow/envs/mapper.yaml`
  - Removed: `workflow/envs/bwa.yaml`, `workflow/envs/bwa_mem2.yaml`

## Config Changes
`config/config_test.yaml` gained the following keys (defaults shown):

```yaml
mapping:
  mapper: bwa         # options: bwa, bwa-mem2
  map_unpaired: true  # include unpaired reads in final BAM
  extra: ""           # extra args to the mapper

outgroups: "data/outgroups.tsv"
outgroups_seq_column: "sequencer"
longread_keywords: ["PacBio", "ONT", "Nanopore"]

# QC and trimming
qc:
  fastqc:
    contaminants: null          # e.g., data/adapters/contaminants.txt
  trimmomatic:
    adapters_fa: null           # e.g., data/adapters/TruSeq3-PE.fa
    clip: "2:30:10"
    slidingwindow: "4:20"
    leading: 3
    trailing: 3
    minlen: 36

# Long-read options
longread:
  filter:
    enabled: true
    min_length: 1000
    keep_percent: 90
    extra: ""
  mapping:
    extra: ""     # extra args passed to minimap2
```

## How It Works
- Short-read outgroups (inferred from `outgroups.tsv:sequencer`):
  - Pre-trim QC: `fastqc_outgroup_pre` → `multiqc_outgroup_pre`
  - Trimming: `trimmomatic_outgroup_pe` (paired + unpaired outputs)
  - Post-trim QC: `fastqc_outgroup_post` → `multiqc_outgroup_post`
  - Mapping: `map_outgroup_paired` (+ `map_outgroup_unpaired` if enabled) → `outgroup_final_bam`
- Long-read outgroups:
  - Merge single-end fastq: `merge_fastq_long` (optional, if input is split per SRR)
  - Optional filtering: `filtlong_outgroup`
  - Mapping: `map_outgroup_long_minimap2` with preset inferred per-sample
- Aggregate target: `outgroup_bams` includes both short- and long-read outgroups.

## Running Examples
- Run short-read QC and trimming for all outgroup samples:
  - `snakemake -s workflow/Snakefile -c 8 outgroup_qc_trim_all`
- Produce final BAMs for all outgroups (short and long):
  - `snakemake -s workflow/Snakefile -c 12 outgroup_bams`

## Notes
- SRR fetch remains in `workflow/rules/outgroup.smk` (`fetch_sra`). Use it locally as per your environment policy. The cluster flow assumes merged fastqs already exist.
- The old example placeholders in `_archive/` or `workflow/rules/fastp.smk`, `workflow/rules/bwa.smk`, etc., are not used by the new flow.
- If you want long-read stats aggregated (e.g., `samtools flagstat` → MultiQC), that can be added easily.
- If desired, `outgroup_bams` can be wired into `rule all` behind a config flag.

