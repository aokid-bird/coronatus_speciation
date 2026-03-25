# Refactoring Steps Completed

This document summarizes the refactoring that has already been applied to the pipeline, in the order it was carried out. It is meant to explain the current structure of the workflow rather than preserve the prompt-by-prompt conversation history.

## 1. Split local-only and cluster-only workflow modules

The first structural change was to stop mixing local bootstrap work and cluster analysis work in the same rule files.

What changed:
- Reference acquisition was split into:
  - `workflow/rules/reference_local.smk`
  - `workflow/rules/reference_shared.smk`
- Outgroup acquisition and shared outgroup metadata were split into:
  - `workflow/rules/outgroup_local.smk`
  - `workflow/rules/outgroup_shared.smk`
- The main `workflow/Snakefile` now includes local-only modules only when `environment: local`, and cluster-only modules only when `environment: cluster`.

Why:
- A rule file should not mix divergent execution environments.
- It is easier to read `rule all` and easier to reason about what is expected on local versus on cluster.

## 2. Reorganize the config templates around grouped workflow responsibilities

The next step was to update the config templates so they reflect the current workflow structure instead of older flat or partially duplicated settings.

What changed:
- `templates/config/defaults_local.yaml`
- `templates/config/defaults_cluster.yaml`

The templates now group settings into sections for:
- project/basic settings
- storage paths
- metadata and sample/outgroup inputs
- reference settings
- QC and mapping
- analysis blocks
- ANGSD arguments
- minInd ratios

Why:
- The config is now easier to scan and easier to keep aligned with the rule grouping.
- The same logical setting is less likely to appear in multiple places.

## 3. Add per-analysis `enabled` switches

Each major analysis block now has an `enabled` key in config.

What changed:
- Shared analysis activation logic was centralized in:
  - `workflow/rules/common_analysis.smk`
- `workflow/Snakefile` now conditionally includes optional analysis modules and target groups based on these flags.

Examples:
- `pcangsd.enabled`
- `ngsadmix.enabled`
- `ngsdist.enabled`
- `raxml.enabled`
- `treemix.enabled`
- `abbababa2.enabled`
- `snapp.enabled`

Why:
- Users can select which downstream analyses to run from config instead of editing workflow code.

## 4. Unify outgroup inclusion and sample selection logic

Several analyses had different ways of selecting ingroups, excluding samples, and including outgroups.

What changed:
- Shared sample and outgroup helper logic was centralized in:
  - `workflow/rules/common_samples.smk`
  - `workflow/rules/common_analysis.smk`
- Analysis-specific bamlist generation now uses shared selection helpers.
- `ngsdist` gained its own independent `include_outgroups` handling instead of inheriting `angsd_global` behavior accidentally.

Why:
- The pipeline now uses one consistent model for:
  - ingroup selection
  - outgroup inclusion
  - sample exclusion
  - sliced outgroup BAM use

## 5. Clean and simplify `rule all` and the main Snakefile

The old top-level workflow had grown difficult to follow.

What changed:
- `workflow/Snakefile` now includes shared modules first, then environment-specific modules.
- Target lists are exported from the rule modules themselves.
- `rule all` is built from grouped target lists rather than hardcoded output paths scattered in one file.

Why:
- The top-level workflow now mostly answers two questions:
  - which modules are active
  - which grouped targets belong in this run

## 6. Reduce `common.smk` into narrower shared modules

The old `common.smk` had accumulated config parsing, metadata loading, path logic, helper functions, and analysis state in one place.

What changed:
- The shared layer was split into:
  - `workflow/rules/common_config.smk`
  - `workflow/rules/common_paths.smk`
  - `workflow/rules/common_samples.smk`
  - `workflow/rules/common_analysis.smk`
  - `workflow/rules/common_rules.smk`
- `workflow/rules/common.smk` is now just a compatibility shim.

Why:
- Shared responsibilities are now separated by role:
  - config normalization
  - filesystem/path setup
  - sample/outgroup metadata
  - analysis activation and shared analysis state
  - generic shared rules

## 7. Standardize target ownership inside rule modules

Several outputs used to be assembled centrally in the Snakefile even though they really belonged to specific modules.

What changed:
- Major rule files now export their own target lists.
- This was done for preprocessing, mapping, and downstream analyses including:
  - `pcangsd.smk`
  - `ngsadmix.smk`
  - `ngsdist.smk`
  - `angsd_raxml.smk`
  - `angsd_treemix.smk`
  - `angsd_sfs.smk`
  - QC and mapping modules

Why:
- The file that produces outputs should define the corresponding target list when possible.
- This keeps the Snakefile thinner and lowers coupling.

## 8. Standardize figures under `figures/`

Figure outputs were split between `results/` and `figures/`, and some plots used PNG while others used PDF.

What changed:
- Main plotted outputs were normalized toward `figures/...`
- Several R script wrappers were adjusted so figure directories are created explicitly.
- A number of analysis plots now follow the RAxML-style figure placement more consistently.

Why:
- Figures are easier to find when they live under one root rather than being mixed into analysis result directories.

## 9. Pull resources and threads into config

Thread and resource settings had been embedded directly inside many rules.

What changed:
- Shared resource/thread resolution now lives in:
  - `workflow/rules/common_config.smk`
- The templates now define top-level sections:
  - `resources`
  - `threads`
- Rule files read these shared settings instead of carrying scattered embedded defaults where practical.

The current layout includes sections such as:
- `resources.reference`
- `resources.qc`
- `resources.mapping`
- `resources.analyses`
- `threads.reference`
- `threads.qc`
- `threads.mapping`
- `threads.analyses`

Why:
- Runtime tuning is now easier to manage from config.
- Rule files remain readable because the config resolution is shared.

## 10. Add support for memory strings like `4G`

The workflow previously expected numeric `mem_mb` values.

What changed:
- `workflow/rules/common_config.smk` now parses memory strings such as:
  - `4G`
  - `16G`
  - `200G`
- These are normalized to `mem_mb` integers before Snakemake uses them.

Why:
- Human-readable config values are easier to maintain, while Snakemake still receives the integer MB values it expects.

## 11. Standardize conda environments for generic R and Python scripts

Several small R and Python helper envs had drifted apart even though they supported the same style of scripts.

What changed:
- Added:
  - `workflow/envs/python_utils.yaml`
  - `workflow/envs/r_plotting.yaml`
- Repointed generic helper-script rules to those envs.
- Removed redundant generic envs:
  - `workflow/envs/python_pandas.yaml`
  - `workflow/envs/intersect_sites.yaml`
  - `workflow/envs/ngsrelate_input.yaml`
  - `workflow/envs/pca_selection.yaml`
  - `workflow/envs/plot.yaml`
  - `workflow/envs/plot_kinship.yaml`
  - `workflow/envs/ngsdist_post.yaml`

What stayed separate:
- `treemix_eval.yaml`
- `plot_tree.yaml`
- `vcfR.yaml`
- `abbababa2_estavg.yaml`
- `snapp_prep.yaml`
- `beast_snapp.yaml`

Why:
- Generic helper and plotting scripts now share standardized environments.
- Specialized tools still keep dedicated envs where dependency surfaces are materially different.

## 12. Standardize rule file headers and improve rule docstrings

Rule files now include clearer file-level headers, especially for the main analysis and preprocessing modules.

What changed:
- Main `.smk` modules now explain:
  - what the file does
  - which inputs it expects
  - which config keys it reads
- Rule-level docstrings were cleaned up where they were too terse or inconsistent.

Why:
- This makes the workflow much easier to navigate for future edits and debugging.

## 13. Validation practice used during this refactor

After each substantial refactoring step, the workflow was revalidated by activating the Snakemake conda environment and running dry-runs with validation configs.

Validation method:
1. Activate the environment described in `docs/codex/run_pythong.sh`
2. Run Snakemake dry-runs against:
   - `config/validation_local.yaml`
   - `config/validation_cluster.yaml`
   - `config/validation_cluster_minimal.yaml`

Note:
- In this Codex environment, `HOME` was redirected to a workspace-local directory during validation because Snakemake tried to write its cache under a location blocked by the sandbox.

## Current workflow shape

At this point, the pipeline has:
- environment-pure rule modules
- grouped config templates
- per-analysis enable switches
- centralized sample/outgroup logic
- centralized resource/thread config
- standardized generic conda envs
- clearer file headers and rule docstrings

The remaining work is mostly incremental cleanup rather than structural redesign.
