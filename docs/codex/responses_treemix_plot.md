# Codex – TreeMix Plot Integration Responses

This document records my responses and actions taken to integrate plotting for TreeMix into the Snakemake pipeline, based on the prompt in `docs/codex/prompts_treemix_plot.md`.

## Prompt Summary

- Integrate plotting and evaluation from an old R script (`codes_from_old_projects/figures_2_rlrmv_mindepthind2_treemix.R`) into the Snakemake pipeline.
- Generalize project-specific variables:
  - Outgroups from config rather than fixed (e.g., Laja).
  - Bamlist from pipeline outputs (not hardcoded).
  - Number of replicates and migration edges from Snakemake config.
  - Use khroma color package; “bright” palette in original order.
- Avoid modifying plotting logic unless necessary; ask if major changes are required.
- Ensure plots save correctly on clusters (avoid ggsave issues).

## Implemented Changes

### New Plotting Script

- Added: `workflow/scripts/treemix_make_plots.R`
  - Reads TreeMix evaluation outputs (`eval_runs.csv`, `eval_summary.csv`).
  - Selects best replicate per migration edge by log-likelihood.
  - Generates three PDF figures under `results/treemix/{output_prefix}/{mode}/plots/`:
    - `treemix_validation2.pdf`: stacked lnL and Variance Explained.
    - `treemix_resid_{mode}.pdf`: residual heatmaps for edges 0..5 with consistent scaling.
    - `treemix_trees_{mode}.pdf`: multi-panel trees (0..5 migration edges) using `ggplot_treemix.R`.
  - Sources `workflow/scripts/treemix_plotting_funcs.R` and `workflow/scripts/ggplot_treemix.R`.
  - Creates a `popordcol` file from `.cov.gz` headers to drive residual/heatmap ordering.
  - Color palette: khroma “bright” in original order.
  - Uses explicit PDF devices for cluster environments; avoids `ggsave`.
  - Generalizes `max_m` and `reps` via arguments from Snakemake/config.

### Snakemake Rule

- Updated: `workflow/rules/angsd_treemix.smk`
  - New rule `treemix_plots` runs the plotting script after evaluation.
  - Outputs (PDF only):
    - `results/treemix/{output_prefix}/{TREEMIX_MODE}/plots/treemix_validation2.pdf`
    - `results/treemix/{output_prefix}/{TREEMIX_MODE}/plots/treemix_resid_{TREEMIX_MODE}.pdf`
    - `results/treemix/{output_prefix}/{TREEMIX_MODE}/plots/treemix_trees_{TREEMIX_MODE}.pdf`
  - Passes `--mode`, `--max_m`, `--reps`, samples TSV, bamlist, plotting functions, and `ggplot_treemix.R`.
  - Adds `KHROMA_LIB` env passthrough from `config.treemix.khroma_lib` for optional manual `khroma` install.

- Updated: `workflow/Snakefile` default targets to include the three PDF plots for TreeMix.

### Environment

- Updated: `workflow/envs/treemix_eval.yaml`
  - Added `r-ggplot2>=3.5.0` (for `ggnewscale`).
  - Added `r-khroma`, `r-ggnewscale`, `r-patchwork`, and `r-cairo`.

### khroma Loading Without Affecting ggplot2

- Problem: Prepending a user library via `.libPaths()` can load older `ggplot2`, causing version warnings.
- Solution in `treemix_make_plots.R`:
  - Do not modify `.libPaths()` globally.
  - Try `library(khroma)` from the conda env first.
  - If not present, load only `khroma` from a specific path via `KHROMA_LIB` environment variable using `library(khroma, lib.loc=...)`.
  - Wrap `khroma::colour()` access in a helper `khroma_colour()` to avoid repeating namespace lookups.
  - Rule `treemix_plots` forwards `config.treemix.khroma_lib` into `KHROMA_LIB`.

## File References

- `workflow/scripts/treemix_make_plots.R:1`
- `workflow/rules/angsd_treemix.smk:1`
- `workflow/Snakefile:1`
- `workflow/envs/treemix_eval.yaml:1`
- `workflow/scripts/treemix_plotting_funcs.R:1`
- `workflow/scripts/ggplot_treemix.R:1`

## Usage Notes

- Ensure Snakemake uses conda: `--use-conda` (and `--conda-frontend mamba` recommended).
- If `khroma` isn’t available in conda on your platform, set a manual library path:
  - In config: `treemix.khroma_lib: "/home/USER/R/x86_64-redhat-linux-gnu-library/4.1"`
  - Or per run: `KHROMA_LIB=/path/to/Rlib snakemake --use-conda ...`
- The plotting rule reads evaluation outputs, so it runs after `treemix_eval`.

## Rationale and Alignment With Prompt

- Generalized outgroups, replicates, and edges via existing pipeline config and rules.
- Embedded plotting logic closely mirroring the previous scripts; only minimal, necessary changes (device handling, color loading, parameterization) were made.
- Saved plots in a cluster-safe way (PDF devices) per your request to avoid PNG/EPS.

## Potential Follow-ups

- If you want additional formats (PNG for quick previews), we can add optional conversions from the PDFs.
- If further customization of the color palette order is needed for specific taxa, we can extend the mapping while keeping khroma as the base palette.

