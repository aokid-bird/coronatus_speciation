# ANGSD → TreeMix Integration (glactools binary, outgroups, LD block size, evaluation)

This page documents how the legacy R script `angsd_treemix_mindepthind2_minind0.5.R` and `lddecay_blocksizefinder.R` were integrated into the Snakemake pipeline. It captures your prompts, Codex responses, and the concrete changes made.

## Your Prompts (condensed)

- General
  - Adapt legacy R code to the current Snakemake pipeline (config, TSVs, dirs).
  - Ignore legacy PBS `cluster.script`; we run on PBS via Snakemake unless `environment == local`.
  - May add columns to existing TSVs for metadata.
  - Replace R loops with Snakemake rules/expands.

- Script 1: `angsd_treemix_mindepthind2_minind0.5.R`
  - 2) ANGSD to BCF; 3) convert BCF→ACF with glactools and specify populations; produce TreeMix input.
  - 3–5) Compute block size parameter for TreeMix via `lddecay_blocksizefinder.R`.
  - 4) Run TreeMix; 5) Evaluate models using `plotting_funcs.R` (TreeMix), OptM and custom functions.
  - Run multiple subsets: ignore; use same pipeline with different configs instead.
  - A new `angsd_treemix.smk` exists: modify it to integrate pipeline.
  - Support outgroups (defined in `outgroups.tsv`) and let users choose: include outgroups or not; per-population vs per-individual analyses.
  - Use the binary glactools (hard to compile) placed at `workflow/bin/glactools`.
  - Should the pipeline include `~/treemix-1.13/src/plotting_funcs.R`? (optional).

- Additional prompts during iteration
  - Outgroup mapping integration exists (short and long reads); use outgroup BAMs.
  - Include TreeMix targets in `rule all`.
  - Handle multiple outgroup BAMs. If `root_label` not set, derive from `taxon` in `outgroups.tsv` by concatenating the first two letters of each word (e.g., “Gallinago media” → “Game”).
  - Prefer to retain separate outgroup taxa (do not merge into a single ‘OUT’), but allow a config switch.

## Codex Response Summary

- Implemented a full TreeMix workflow integrated with ANGSD and ngsLD outputs.
- Added binary glactools usage (no compilation needed) via `workflow/bin/glactools` (path configurable).
- Added robust outgroup handling including per-population and per-individual modes, outgroup inclusion toggles, and root selection rules driven by `outgroups.tsv`.
- Kept evaluation modular: OptM optional; TreeMix `plotting_funcs.R` optional via config.

## Files Added or Updated

- New/Updated Rules
  - `workflow/rules/angsd_treemix.smk`
    - `make_bamlist_treemix`: combines unrelated ingroup BAMs with outgroup BAMs (optional) into a treemix bamlist.
    - `angsd_treemix`: ANGSD on intersecting sites to BCF/geno.
    - `bcf_to_acf`: bcftools convert/annotate (drop FORMAT/GL) → glactools `vcfm2acf` and `index` using the binary.
    - `glactools_prepare_treemix`: Python helper to rename, merge populations, set root (outgroups), and write `input.treemix` from ACF.
    - `treemix_input_gz`: compresses input.treemix.
    - `lddecay_blocksize`: runs R block-size finder (standalone), writes `lddecay/block_size.txt` + plot.
    - `treemix_run`: runs TreeMix for edges 0..`treemix.max_m` and replicates 1..`treemix.reps` with `-k` from block size and `-root` when appropriate.
    - `treemix_eval`: summarizes lnL across runs; optionally uses OptM and TreeMix `plotting_funcs.R`.

- New Scripts
  - `workflow/scripts/prepare_glactools_treemix.py`
    - Inputs: `--acf`, `--bamlist`, `--samples` (TSV with `sample` and group column), `--group-col`, `--mode pop|indiv`.
    - Outgroups: `--include-outgroups`, `--outgroups-tsv`, `--root-label` (optional), `--merge-outgroups true|false`.
    - Uses glactools: `view -p` (best-effort sample order check), `rename`, `meld`, `usepopsrootanc`, `replaceanc`, `acf2treemix`.
    - Outgroup root inference:
      - If `root_label` in config: use it.
      - Else if merging outgroups: single label based on `taxon` if unique else `OUT`.
      - Else (retain separate taxa): derive labels from `taxon` (e.g., “Game”, “Gaga”) and pick the taxon with most present samples as root.
  - `workflow/scripts/lddecay_blocksizefinder.R`
    - Standalone R script to compute block size with LOESS-smoothed LD decay using ACF and ngsLD outputs. Outputs CSVs and `block_size.txt`.
  - `workflow/scripts/treemix_optm_eval.R`
    - Summarizes `.llik` files across edges/replicates, writes `eval_runs.csv` and `eval_summary.csv`.
    - Optionally loads TreeMix `plotting_funcs.R` (for variance explained) and runs OptM if available.

- Config Updates
  - `config/config_test.yaml`
    - Added `treemix` section:
      - `include_outgroups: true|false`
      - `mode: pop|indiv`
      - `glactools_bin: workflow/bin/glactools`
      - `root_label: <string or null>`
      - `max_m: 6`, `reps: 10`
      - `use_optm: true|false`, `plotting_funcs: <path or null>`
      - `merge_outgroups: false` (retain separate outgroup taxa; set true to merge all outgroups into one root group)

- Common Variables
  - `workflow/rules/common.smk`
    - Added `TREEMIX_*` and `GLACTOOLS_BIN` variables derived from config.

- Pipeline Wiring
  - `workflow/Snakefile`
    - Included `rules/angsd_treemix.smk`.
    - Added TreeMix outputs to `rule all`:
      - `results/treemix/{output_prefix}/{TREEMIX_MODE}/eval_summary.csv`
      - `results/treemix/{output_prefix}/{TREEMIX_MODE}/eval_runs.csv`

- Environments
  - `workflow/envs/python_pandas.yaml` (for the Python prep script).
  - `workflow/envs/treemix_eval.yaml` (R + OptM for evaluation; LD script also runs fine in this env).

## Glactools Binary Usage

- Path defaults to `workflow/bin/glactools` (ignored by Git per `.gitignore`).
- Override via `treemix.glactools_bin` in config.
- `lddecay_blocksize` rule augments PATH (`workflow/bin`) so `glactools view` works inside the R script.

## Outgroup Behavior

- Controlled by config:
  - `treemix.include_outgroups`: include outgroup BAMs in bamlist and subsequent steps.
  - `treemix.mode`:
    - `pop`: ingroup merged by `config.group_col`; outgroups merged by taxon (or single label when merging).
    - `indiv`: ingroup remain individuals; outgroup individuals merged by taxon (or into one group when merging).
  - `treemix.merge_outgroups`:
    - `false` (default here): retain separate taxa (e.g., Game, Gaga). Root is chosen or provided.
    - `true`: collapse outgroups into one root population (`treemix.root_label` if set; else derived/“OUT”).
  - `treemix.root_label`: optional explicit root group name; otherwise derived as above.

## Usage

- Ensure outgroup BAMs exist (from mapping_outgroup.smk flow).
- Place glactools binary at `workflow/bin/glactools` (chmod +x), or set `treemix.glactools_bin`.
- Optionally set `treemix.plotting_funcs` (path to TreeMix `plotting_funcs.R`) and `treemix.use_optm`.
- Run full pipeline (including TreeMix):
  - `snakemake -s workflow/Snakefile -c 12 all`
- Or just TreeMix evaluation targets:
  - `snakemake -s workflow/Snakefile -c 12 results/treemix/{output_prefix}/{TREEMIX_MODE}/eval_summary.csv`

## Outputs

- `results/treemix/{output_prefix}/{mode}/input.acf.gz`
- `results/treemix/{output_prefix}/{mode}/input.treemix[.gz]`
- `results/treemix/{output_prefix}/{mode}/lddecay/block_size.txt` (+ CSV + plot)
- `results/treemix/{output_prefix}/{mode}/treemix_e{edge}_o{rep}.{llik,vertices.gz,edges.gz}`
- `results/treemix/{output_prefix}/{mode}/eval_runs.csv`, `eval_summary.csv`, optional `eval_optm.tsv`/`eval_optm.pdf`

## Notes

- The ACF is rooted via glactools before TreeMix; `treemix_run` passes `-root` only if `treemix.root_label` is explicitly set.
- If `outgroups.tsv` lacks `taxon`, separate-taxa mode falls back to per-sample labels for outgroups.
- The LD half-life approach is a practical approximation (LOESS) mirroring your logic but without external R sources; can be swapped for a parametric fit if desired.

## Quick Change Log

- New rules: TreeMix pipeline (`angsd_treemix.smk`).
- New scripts: `prepare_glactools_treemix.py`, `lddecay_blocksizefinder.R`, `treemix_optm_eval.R`.
- Config: new `treemix` section with binary path, toggles, and parameters.
- Rule all: TreeMix evaluation added.
- Envs: `python_pandas`, `treemix_eval`.

## Addendum: Sliced Outgroups Integration (This Session)

- Switched TreeMix bamlist to use sliced outgroup BAMs at `{bam_dir}/outgroups_sliced/` generated from the union of ingroup sites, reducing IO and ensuring site compatibility with intersect filters.
- Added dynamic `-minInd` to `angsd_treemix` computed from the total number of individuals in the actual bamlist (including outgroups if enabled), based on `minIndRatio.treemix` (fallback to `global`/`intersect`).
- Related cross-pipeline updates in this session:
  - Union site builder and outgroup BAM slicing rules were added (`workflow/rules/slice_outgroups.smk`, `workflow/scripts/union_sites.py`).
  - Per-analysis toggles for outgroup inclusion added (`angsd_global.include_outgroups`, `angsd_raxml.include_outgroups`; Treemix keeps `treemix.include_outgroups`).
  - Plotting metadata now includes outgroup taxa as separate populations via `taxon` in `outgroups.tsv`.

## Full Transcript (verbatim excerpts)

User Prompt:
"I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder \"codes_from_old_projects\", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? ... 1) angsd_treemix_mindepthind2_minind0.5.R ... I have prepared a smk file called \"angsd_treemix.smk\". Please modify this for pipeline adaptation. ... I want to include ougroups to this analysis. ... Please use lddecay_blocksizefinder.R ... \"~/treemix-1.13/src/plotting_funcs.R\" is provided by treemix. Should the pipeline include this?"

Codex Response:
"I’m going to scan the repo for the treemix and ANGSD pieces, then open the related R scripts to design Snakemake rules and config hooks."

System actions: Searched repo, opened `angsd_treemix.smk`, `angsd_global.smk`, `common.smk`, R scripts, mapping rules.

User Prompt:
"Binary version of the Glactools will be used since compilation of the package is difficult in the environment."

Codex Response:
"I’m going to add treemix config defaults, create helper scripts, and expand angsd_treemix rules to use the glactools binary, outgroup toggles, and block-size + treemix evaluation."

System actions (edits):
- Added TREEMIX config vars to `workflow/rules/common.smk`.
- Added `workflow/scripts/prepare_glactools_treemix.py`.
- Added `workflow/scripts/lddecay_blocksizefinder.R` (standalone).
- Added `workflow/scripts/treemix_optm_eval.R`.
- Expanded `workflow/rules/angsd_treemix.smk` to full flow (ANGSD→VCF→ACF→prep→blocksize→runs→eval) using `workflow/bin/glactools`.
- Updated `config/config_test.yaml` with `treemix` section.
- Included `rules/angsd_treemix.smk` in `workflow/Snakefile`.
- Added `workflow/envs/treemix_eval.yaml`.

Codex Summary to User:
- Described integrated treemix flow, glactools binary usage, outgroup handling, scripts, and how to run. Asked if results should be wired into rule all.

User Prompt:
"Yes include it in the rule all target. Also, I was uncertain how it works when there are more than two outgroup bams. And in case root_label not specified in the config, please use the \"taxon\" column of the outgroups.tsv where two first letters of each word are glued for each taxon (e.g. Gallinago media -> Game, Gallinago gallinago -> Gaga)."

Codex Response:
"I’ll wire TreeMix targets into rule all and enhance outgroup-root handling as requested."

System actions (edits):
- `workflow/Snakefile`: added TreeMix eval targets to `rule all`.
- `workflow/rules/common.smk`: allowed `root_label` to be null.
- `prepare_glactools_treemix.py`: derive outgroup labels from `taxon` when root_label is not provided; support multiple outgroups; ensure ACF rooted; indiv/pop modes preserve outgroup logic.
- Added `workflow/envs/python_pandas.yaml` and switched the prep rule to use it.

User Prompt:
"If you want TreeMix to retain separate outgroup taxa (instead of merging to “OUT”), we can add a config switch ... -> Yes, I want to retain them separate."

Codex Response:
"Great — I added a switch to retain separate outgroup taxa and wired it through the pipeline."

System actions (edits):
- `config/config_test.yaml`: added `treemix.merge_outgroups: false`.
- `prepare_glactools_treemix.py`: added `--merge-outgroups` and implemented per-taxon vs merged behavior.
- `angsd_treemix.smk`: passed `--merge-outgroups` from config.

User Prompt:
"Please make a document called \"docs/codex/treemix_angsd.md\" that include 1) my prompts, 2) your responses and descriptions on changes, and everything you have done in this chat. Then please update docs/codex/index.md."

Codex Response:
"Added comprehensive TreeMix+ANGSD docs and wired them into the index. If you want the full raw prompt/response transcripts instead of condensed summaries, I can append them verbatim."

User Prompt:
"Yes, Please do so"

Codex Response:
"Appended this transcript section with verbatim prompt/response excerpts and edit summaries."
