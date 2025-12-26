# Including Sliced Outgroups Across ANGSD and Downstream Analyses

This document records the prompt(s), the responses, and the concrete changes made to the Snakemake pipeline to support including outgroups via sliced BAMs, while keeping downstream analyses restricted to ingroup-derived sites. It excludes any response about not being able to conduct a dry-run in the current sandbox.

## Prompts
- Develop rules to slice raw outgroup BAMs to the union of sites obtained from ingroups (full join A ∪ B ∪ C over per‑population angsd_intersect outputs), to make ANGSD on outgroups lighter.
- Keep `angsd_intersect.smk` ingroup‑only; outgroups are not part of intersect site discovery.
- Let users decide whether to include outgroups separately for: 1) `angsd_global.smk`, 2) `angsd_raxml.smk`, 3) `angsd_treemix.smk`. When included, use sliced outgroup BAMs. Ensure filters like `-minInd` reflect the total number of individuals given the inclusion/exclusion choice.
- Update `ngsadmix.smk` and `pcangsd.smk` so outgroup taxa are treated as distinct populations named by `taxon` in `data/outgroups.tsv`.

## What I Implemented
- Added a slicing workflow that constructs the union of ingroup sites and slices outgroup BAMs to those sites.
- Kept ingroup‑only intersect logic unchanged; downstream ANGSD rules consume ingroup intersect sites, optionally with sliced outgroup BAMs included in the bamlist.
- Added per‑analysis toggles to include outgroups for global, RAxML, and TreeMix flows. `-minInd` is computed based on the actual bamlist used in each analysis.
- Augmented metadata used for plotting (PCA, Admixture) so outgroups appear as separate populations using their `taxon` labels.

## Detailed Changes
- New union + slicing rules
  - `workflow/rules/slice_outgroups.smk`
    - `union_sites_groups`: Builds A ∪ B ∪ C from `results/angsd_intersect/{output_prefix}/{group}/gl.geno.gz` and writes:
      - `results/union_sites/{output_prefix}/union.sites` (TSV: scaf pos)
      - `results/union_sites/{output_prefix}/union.chr`
      - `results/union_sites/{output_prefix}/union.bed` (0‑based intervals for `samtools view -L`)
    - `slice_outgroup_bam`: Slices `{bam_dir}/outgroups/{sample}.bam` to union `.bed` and outputs `{bam_dir}/outgroups_sliced/{sample}.bam` + `.bai`.
    - `outgroup_bams_sliced`: Aggregate target to ensure sliced BAMs exist for all outgroup samples.
  - `workflow/scripts/union_sites.py`: Helper to create union site list and scaffolds from multiple ANGSD `.geno.gz` files.
- Common parameters and toggles
  - `workflow/rules/common.smk`
    - `OUTGROUP_SLICED_DIR = {bam_dir}/outgroups_sliced`.
    - Per‑analysis toggles:
      - `ANGSD_GLOBAL_INCLUDE_OUTGROUPS` from `angsd_global.include_outgroups` (bool; default false).
      - `ANGSD_RAXML_INCLUDE_OUTGROUPS` from `angsd_raxml.include_outgroups` (bool; default false).
    - Helper `get_minInd_ratio(key)` with fallbacks to `minIndRatio.global` or `minIndRatio.intersect`.
- ANGSD global + LD + unrelated/unlinked
  - `workflow/rules/angsd_global.smk`
    - `make_bamlist_global_analysis`: starts from ingroup `global` bamlist; appends sliced outgroups when enabled.
    - `angsd_global`: uses analysis bamlist; appends `-minInd` = floor(total_samples × `minIndRatio.global`).
    - `ngsld_global`: uses same analysis bamlist for consistent `n_ind`.
    - `make_bamlist_unrelated_analysis`: starts from unrelated ingroup; appends sliced outgroups when enabled.
    - `angsd_global_unrelated_unlinked`: uses unrelated analysis bamlist; appends `-minInd` similarly.
- RAxML flow
  - `workflow/rules/angsd_raxml.smk`
    - `make_bamlist_raxml_analysis`: starts from unrelated ingroup; appends sliced outgroups when enabled.
    - `angsd_raxml`: uses that bamlist; appends `-minInd` from `minIndRatio.raxml` (fallbacks to global/intersect).
    - Downstream `catg_format`, `plot_raxml` now consume this bamlist.
- TreeMix flow
  - `workflow/rules/angsd_treemix.smk`
    - `make_bamlist_treemix`: now appends sliced outgroup BAMs (respects `treemix.include_outgroups`).
    - `angsd_treemix`: appends `-minInd` from `minIndRatio.treemix` (fallbacks to global/intersect).
- PCAngsd & NGSadmix plotting metadata
  - `workflow/rules/ngsadmix.smk` and `workflow/rules/pcangsd.smk`
    - New helper rules create `results/metadata/{output_prefix}/samples_plus_outgroups.tsv` by appending rows for outgroup `sample_id` with `group_col` set from their `taxon`.
    - Plotting rules now read the augmented metadata so outgroups appear as separate populations.
- Snakefile wiring
  - `workflow/Snakefile`: includes `rules/slice_outgroups.smk` and adds union/sliced outputs to `rule all` (cluster mode).

## Config Additions and Defaults
- `minIndRatio` keys added in example config: `global`, `raxml`, `treemix` (defaults 0.5 in `config/config_test.yaml`).
- Per‑analysis inclusion toggles:
  - `angsd_global.include_outgroups: false` (default)
  - `angsd_raxml.include_outgroups: false` (default)
  - `treemix.include_outgroups: true|false` (existing; unchanged)

## Usage
- Intersect sites remain ingroup‑only via `angsd_intersect.smk`.
- Build sliced outgroup BAMs based on union of ingroup sites:
  - Generated automatically when downstream analyses depend on `{bam_dir}/outgroups_sliced/{sample}.bam` and union site outputs.
- Enable outgroups per analysis by config:
  - Global ANGSD: set `angsd_global.include_outgroups: true`.
  - RAxML: set `angsd_raxml.include_outgroups: true`.
  - TreeMix: set `treemix.include_outgroups: true`.
- `-minInd` is computed from the analysis bamlist size using `minIndRatio.<analysis>`. If the key is missing, falls back to `minIndRatio.global` or `minIndRatio.intersect`.

## Notes
- Sliced outgroup BAMs are produced at `{bam_dir}/outgroups_sliced/`.
- Plotting color order still uses `config.populations`; to control legend/ordering for outgroup taxa, add their labels (matching `taxon` abbreviations or names) to that list if desired.
- `angsd_common_args` and `angsd_args[...]` are unchanged; `-minInd` is appended by rules when applicable.

## Files Touched
- Added
  - `workflow/rules/slice_outgroups.smk`
  - `workflow/scripts/union_sites.py`
- Updated
  - `workflow/rules/common.smk`
  - `workflow/rules/angsd_global.smk`
  - `workflow/rules/angsd_raxml.smk`
  - `workflow/rules/angsd_treemix.smk`
  - `workflow/rules/ngsadmix.smk`
  - `workflow/rules/pcangsd.smk`
  - `workflow/Snakefile`
  - `config/config_test.yaml`

## Future Improvements
- Optionally cache union sites per scenario to avoid recomputation when only outgroup choices change.
- Expand plotting to automatically include outgroup taxa in `config.populations` or derive palettes from `samples_plus_outgroups.tsv`.
- Add validation rules to ensure outgroup BAMs exist (local vs cluster workflows) and union slicing completes before ANGSD starts.

