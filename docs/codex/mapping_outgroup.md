# Codex Integration Notes: mapping_outgroup.smk

This document focuses on the mapping integration for outgroups, including short-read (bwa/bwa-mem2) and long-read (minimap2) flows, with consolidation of mapping environments.

## User Prompts and Codex Responses (mapping-related)

### Prompt (User)
I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake.

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

### Response (Codex)
- Implemented mapping rules for outgroups that:
  - Choose bwa or bwa-mem2 based on config
  - Optionally map unpaired reads and merge with paired
  - Write sorted/indexed BAMs to `config['bam_dir']/outgroups`
- Added long-read branch using minimap2, with optional filtlong filtering
- Derived read type and minimap2 presets from outgroup metadata in `outgroups.tsv`

### Prompt (User)
Please do or think about the following
- bwa.yaml and mapper.yaml seem to be redundant. Please consolidate them and reflect changes in rules that use these environments
- add minimap2-based trimming/mapping for long reads and integrate those outputs as well

### Response (Codex)
- Consolidated mapping envs into one `workflow/envs/mapper.yaml` (bwa, bwa-mem2, samtools, minimap2, filtlong)
- Removed `workflow/envs/bwa.yaml` and `workflow/envs/bwa_mem2.yaml`
- Switched reference indexing rules to use consolidated env
- Added long-read filtering (filtlong) and mapping (minimap2) rules and integrated them into outgroup BAM aggregation

---

## What mapping_outgroup.smk does
File: `workflow/rules/mapping_outgroup.smk`

- Short-read mapping (from trimmed pairs):
  - `map_outgroup_paired`
    - Inputs: `results/trimmomatic/outgroup/{sample}_pair_R1/2.fastq.gz`, reference fasta and mapper index
    - Mapper: `bwa mem` or `bwa-mem2 mem` selected via `config.mapping.mapper`
    - Output: temporary paired BAM in `results/mapping/{output_prefix}/outgroup/{sample}.paired.bam`
  - `map_outgroup_unpaired` (optional, if `config.mapping.map_unpaired: true`)
    - Inputs: `results/trimmomatic/outgroup/{sample}_unpair_R1/2.fastq.gz`
    - Output: temporary unpaired BAMs `results/mapping/{output_prefix}/outgroup/{sample}.unpaired_R{read}.bam`
  - `outgroup_final_bam`
    - Merges paired (+unpaired if enabled), sorts, and indexes to `{config_bam_dir}/outgroups/{sample}.bam`

- Long-read mapping (single-end):
  - `filtlong_outgroup` (optional; controlled by `config.longread.filter.enabled`)
    - Input: `data/merged/outgroup/{sample}.fastq.gz`
    - Output: `results/longread/filter/{sample}.fastq.gz`
  - `map_outgroup_long_minimap2`
    - Input: filtered fastq (or merged fastq if filtering disabled), reference
    - Mapper: `minimap2` with preset inferred per-sample from `sequencer` column
      - `map-ont` for ONT/Nanopore, `map-hifi` if “HiFi”, else `map-pb`
    - Output: `{config_bam_dir}/outgroups/{sample}.bam` + `.bai`

- Aggregate target:
  - `outgroup_bams`: ensures all short- and long-read outgroup BAMs are present

### Supporting directories and variables
- Trimmed pairs: `results/trimmomatic/outgroup`
- Long-read filtered: `results/longread/filter`
- Intermediate mapping: `results/mapping/{output_prefix}/outgroup`
- Final BAM destination: `{config['bam_dir']}/{output_prefix}/outgroups`
- Mapper backend and reference paths come from the shared config/path layer:
  - `MAPPER`, `REF_UNZIPPED`, `REF_MAP_INDEX_MAIN`

## Metadata-driven behavior
- `workflow/rules/common_samples.smk` reads `data/outgroups.tsv` and infers for each `sample_id`:
  - `OUTGROUP_READ_TYPE`: `short` or `long` based on `outgroups_seq_column` and `longread_keywords`
  - `OUTGROUP_MINIMAP2_PRESET`: minimap2 preset per sample
- Short-read rules operate only on samples with `OUTGROUP_READ_TYPE == 'short'`
- Long-read rules operate on `OUTGROUP_READ_TYPE == 'long'`

## Conda env consolidation
- `workflow/envs/mapper.yaml` contains all mapping tools (bwa, bwa-mem2, minimap2, filtlong, samtools)
- Reference indexing switched to use `../envs/mapper.yaml`
- Removed redundant envs `workflow/envs/bwa.yaml` and `workflow/envs/bwa_mem2.yaml`

## Config keys
`config/config_test.yaml` controls mapping behavior:

```yaml
mapping:
  mapper: bwa          # options: bwa, bwa-mem2
  map_unpaired: true   # include unpaired mapping in final BAM merge
  extra: ""            # extra args to bwa/bwa-mem2

# Long-read options
longread:
  filter:
    enabled: true
    min_length: 1000
    keep_percent: 90
    extra: ""
  mapping:
    extra: ""          # extra args passed to minimap2

# Outgroup metadata
outgroups: data/outgroups.tsv
outgroups_seq_column: sequencer
longread_keywords: [PacBio, ONT, Nanopore]
```

## Usage examples
- Short-read mapping to final BAMs for all outgroups (after trimming):
  - `snakemake -s workflow/Snakefile -c 12 outgroup_bams`
- Long-read mapping (minimap2), with filtering on:
  - Ensure `data/merged/outgroup/{sample}.fastq.gz` exists (or use `merge_fastq_long` from `outgroup_local.smk` if needed)
  - `snakemake -s workflow/Snakefile -c 12 outgroup_bams`

## Notes
- `mapping.extra` is appended to bwa/bwa-mem2 commands; `longread.mapping.extra` to minimap2.
- If you want mapping QC summaries, we can add `samtools flagstat` and a MultiQC aggregation step.
- The aggregate `outgroup_bams` does not change `rule all`; it remains opt-in unless requested otherwise.
