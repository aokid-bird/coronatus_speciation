I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake. Please do this conversion only for the files specified below, since some of the scripts are already converted and integrated into this pipeline.

0) General
Descriptions & Specifications:
- These codes were developed in different projects. Please adapt them to the present snakemake pipelines on config parameters, tsv files to be read, directories, etc.
- cluster.script in the old script intends to convert shell commands to PBS cluster. Please ignore this step since snakemake is intended to be ran in PBS cluster in most cases (except when environment == "local").
- They are ran under environment == "cluster"
- You may add columns to any already-present tsv files to retrieve metadata of samples etc.
- Please adjust styles to snakemake. For example, loop inside R is not adequate but single rules with expand comand would be sufficient in this pipeline.

1) ngsdist_neighbornet_mindepthind2.R (= previous script)
- Description
This script runs ngsDist to calculate pairwise genetic distances among samples. This utilizes a previously-run ANGSD output (beagle) and its bamlist as inputs. 
- Specification
1. For ngsDist, Apptainer/Singularity needs to be used to create a container. The source is here https://github.com/fgvieira/ngsDist. Please implement methods to generate its sif file in lima on Mac Apple Silicon environment. Please make sure that this sif will be transferred to Linux Redhat Cluster. Therefore, make sure that lima is based on x86_64 environment.
2. The ngsDist will be based on rule.angsd_global with its bamlist.
3. Please let users to choose what evolutionary models to be used. For other parameters, for now, please make them as the previous script defines. Although the previous script contains runs for p-distance and JC69, please include one evolutionary model for one rule.
4. dist2nexus is defined in codes_from_old_projects/ngs_functions.R. Please use this R function, but plase determine where to store this function by considering the snakemake piepeline structure.
5. The previous script actually saves nexus by dist2nexus, and then SplitsTree (no CLI) will load this file manually to output nnet_{evolutionary_model}.nexus. Although the previous script tries to plot Neighbornet by loading this file and using phangorn, please let the pipeline to stop by outputting the nexus by dist2nexus.  

---

# Conversation Tracker — ngsDist Integration

This section records the prompts (specifications) and the corresponding Codex responses/implementations, with references to files and outputs in this repository.

## Prompt 1 (Initial Specification)

Summary: Integrate ngsDist using ANGSD outputs; run via Apptainer/Singularity; allow choosing evolutionary models; use `dist2nexus` and stop at writing NEXUS files (SplitsTree usage remains manual).

Response (implemented):
- Container build:
  - Definition: `workflow/singularity/ngsdist.def`
  - Mac Apple Silicon guide (Lima x86_64 + Apptainer): `docs/ngsdist_container_mac.md`
- Rule inputs: Based on `rule angsd_global` outputs and bamlist.
- Rules added: `workflow/rules/ngsdist.smk`
  - `ngsdist_prepare_inputs`: builds `labels.txt`, `posinfo.tsv`, `nsites.txt` from `gl.beagle.gz` and bamlist.
  - `ngsdist_run` (parameterized by evolutionary model; see Prompt 2 changes): runs ngsDist and emits distance matrix.
  - `ngsdist_nexus`: converts ngsDist matrix to NEXUS.
- R helper function location:
  - `dist2nexus` placed at `workflow/scripts/ngs_functions.R` and used by `workflow/scripts/ngsdist_to_nexus.R`.
- Configuration knobs:
  - `config/config_test.yaml`: `ngsdist.models: ["p", "jc69"]`
- Pipeline wiring:
  - `workflow/Snakefile`: includes `rules/ngsdist.smk` and adds ngsDist outputs to `rule all` targets conditionally based on configured models.
- Outputs produced:
  - Distance matrices and NEXUS files under `results/ngsdist_global/{output_prefix}/{model}/`.

Notes:
- The pipeline stops at NEXUS generation as requested. SplitsTree usage remains manual.

## Prompt 2 (Additional Request)

Summary: Consolidate separate p-distance and JC69 rules into a single rule; users select models by names in config; place results under `results/ngsdist_global/{output_prefix}/{evol_model}`.

Response (implemented):
- Consolidated rule:
  - `rule ngsdist_run` with wildcard `{model}`; supported names: `p` and `jc69`.
  - Model mapping: `p -> --evol_model 0`, `jc69 -> --evol_model 1`.
  - Uses `--pos` only for JC69.
- Output directories:
  - Matrix: `results/ngsdist_global/{output_prefix}/{model}/ngsdist`
  - NEXUS: `results/ngsdist_global/{output_prefix}/{model}/ngsdist_input.nexus`
- Snakefile targets:
  - Built dynamically from `ngsdist.models` and point to the per-model directories.

Key file references:
- workflow/rules/ngsdist.smk
- workflow/Snakefile
- config/config_test.yaml
- workflow/scripts/ngs_functions.R
- workflow/scripts/ngsdist_to_nexus.R
