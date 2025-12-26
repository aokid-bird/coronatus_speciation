I am developing a snakemake pipeline for genomics biofinformatics and analyses. I have generated rule files called "angsd_raxml.smk" and "angsd_global.smk" that generally conducts ANGSD and runs downstream analyses. However, in these files, inclusion/exclusion of outgroups is not considered. Furthermore, while "angsd_treemix.smk" does consider inclusion/exclusion of outgroups, this was generated to include the raw outgroup bams. However, this pipeline focuses on the use of intersected sites among the ingroup populations under studies (defined in "angsd_intersects.smk"). Therefore, it needs to make sure that sites used in the downstream analyses are above filtering thresholds (e.g., minIndDepth, minInd) even after including outgroups.
In order to include outgroups in these files and use only intersecting sites of the ingroups, please consider the following and update the pipeline.

1) Please make rules to slice the raw bams (through a flow SRA -> QC/Trims -> Mapper -> BAM) to include obtained sites in ingroups. This does not mean to include "intersecting sites" but rather every site obtained from all the populations. (i.e., if you assume three populations A, B, C, and conducted angsd_intersect for each population, then it will be the full join of the three population sites, i.e., AUBUC in the Venn Diagram). The slice rules are to make bam files as light as possible since reading the raw outgroup BAMs are too huge for ANGSD analyses.

2) "angsd_intersect.smk" does not need to include outgroups. This set of rules is ran only among populations of the ingroup.

3) Let users to select if to add outgroups in the downstream analyses of 1: angsd_global.smk, 2: angsd_raxml.smk, 3: angsd_treemix.smk, but separately. And modify these rules to include the sliced outgroup bams from the ANGSD runs. Filtering parameters such as "minInd" define the states when outgroups are included/excluded. For example, if minIndratio in the config file is set "0.5", then 50% of the all the samples in the bamlist have to have the sites (regardless of the inclusion of outgroups) 

4) Please modify "ngsadmix.smk" and "pcangsd.smk" to be adapted to the modification done in the above. Specifically, the outgroup taxa are to be treated as separate populations named as "taxon" in the outgroups.tsv.

Please first make a document called "docs/codex/including_sliced_ougroup.md" that include 1) my prompts, 2) your responses and descriptions on changes, and everything you have done in this chat, except for your response about that you couldnt conduct dry-run. Then please update docs/codex/index.md.

Implemented changes:
- New rules file `workflow/rules/slice_outgroups.smk`:
  - `union_sites_groups`: builds A ∪ B ∪ C union sites from `results/angsd_intersect/{output_prefix}/{group}/gl.geno.gz`, outputs `results/union_sites/{output_prefix}/union.sites`, `.chr`, and `.bed`.
  - `slice_outgroup_bam`: slices each mapped outgroup BAM (`{bam_dir}/outgroups/{sample}.bam`) to the union `.bed`, writing to `{bam_dir}/outgroups_sliced/{sample}.bam` and index.
  - `outgroup_bams_sliced` aggregate target to ensure all are sliced.
- `angsd_global.smk`:
  - `make_bamlist_global_analysis`: builds analysis bamlist and optionally appends sliced outgroups when `angsd_global.include_outgroups: true`.
  - `angsd_global` and `angsd_global_unrelated_unlinked` now read those bamlists and add `-minInd` computed from `minIndRatio.global` (fallback to `intersect` if not set) based on total individuals in the bamlist.
  - Downstream `ngsld_global` uses the same analysis bamlist.
- `angsd_raxml.smk`:
  - `make_bamlist_raxml_analysis`: starts from unrelated ingroup and optionally appends sliced outgroups when `angsd_raxml.include_outgroups: true`.
  - `angsd_raxml` consumes that bamlist and adds a `-minInd` derived from `minIndRatio.raxml` (fallback to `global`/`intersect`).
- `angsd_treemix.smk`:
  - TreeMix bamlist now uses sliced outgroup BAMs and adds optional `-minInd` derived from `minIndRatio.treemix` (fallbacks same as above). Existing `treemix.include_outgroups` is respected.
- `ngsadmix.smk` and `pcangsd.smk`:
  - Added small helper rules to create `results/metadata/{output_prefix}/samples_plus_outgroups.tsv` by appending rows for outgroup `sample_id` with `group_col` set from `outgroups.tsv` `taxon`. Plotting rules now read this augmented metadata so outgroups appear as separate populations.

Config additions (optional but recommended):
- `angsd_global.include_outgroups: true|false` (default false)
- `angsd_raxml.include_outgroups: true|false` (default false)
- `treemix.include_outgroups: true|false` (already present; continues to work)
- `minIndRatio.global`, `minIndRatio.raxml`, `minIndRatio.treemix` can be set to control `-minInd` by total samples in bamlist. If not set, they fall back to `minIndRatio.intersect` when available.

Notes:
- `angsd_intersect.smk` remains unchanged and only uses ingroup populations.
- Sliced outgroup BAMs live in `{bam_dir}/outgroups_sliced/` and are used whenever outgroups are included.
