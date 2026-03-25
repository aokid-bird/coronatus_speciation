# Codex Integration Docs Index

Welcome. These pages summarize how legacy R scripts were integrated into the Snakemake pipeline for outgroup QC, trimming, and mapping.

- Refactoring Summary: docs/codex/refactoring_steps_completed.md
  - Step-by-step description of the structural refactor already applied across rule grouping, config layout, resources/threads, envs, and documentation.

- QC + Trimming Overview: docs/codex/qc_trimmomatic.md
  - FastQC/MultiQC (pre/post), Trimmomatic parameters, config keys, usage.
- Mapping (Short/Long Reads): docs/codex/mapping_outgroup.md
  - BWA/BWA-MEM2 for short reads, optional unpaired merge, minimap2 + filtlong for long reads, env consolidation, config and examples.
- ANGSD → TreeMix: docs/codex/treemix_angsd.md
  - ANGSD intersect → glactools (binary) → TreeMix input, outgroup handling (pop/indiv, merge-or-separate taxa, root selection), LD-based block size, TreeMix runs, OptM/plotting_funcs evaluation, config and usage.
- Including Sliced Outgroups: docs/codex/responses_including_sliced_ougroup.md
  - Build union of ingroup sites, slice outgroup BAMs, include per-analysis (global/RAxML/TreeMix), dynamic minInd, and plotting metadata for outgroup taxa.

Tips
- Edit config at: config/config_test.yaml
- Core rule files:
  - workflow/rules/qc_trimmomatic.smk
  - workflow/rules/mapping_outgroup.smk
  - workflow/rules/common_config.smk
  - workflow/rules/common_paths.smk
  - workflow/rules/common_samples.smk
  - workflow/rules/common_analysis.smk
  - workflow/rules/common_rules.smk
  - workflow/rules/outgroup_local.smk
  - workflow/rules/outgroup_shared.smk
  - workflow/rules/reference_local.smk
  - workflow/rules/reference_shared.smk
- Mapping tools env: workflow/envs/mapper.yaml

If you want these targets included in the default rule, we can wire an opt-in config flag and update Snakefile accordingly.
