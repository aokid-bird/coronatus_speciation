I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake. Please do this conversion only for the files specified below, since some of the scripts are already converted and integrated into this pipeline.

0) General
Descriptions & Specifications:
- These codes were developed in different projects. Please adapt them to the present snakemake pipelines on config parameters, tsv files to be read, directories, population/sample names etc.
- cluster.script in the old script intends to convert shell commands to PBS cluster. Please ignore this step since snakemake is intended to be ran in PBS cluster in most cases (except when environment == "local").
- They are ran under environment == "cluster"
- You may add columns to any already-present tsv files to retrieve metadata of samples etc.
- Please adjust styles to snakemake. For example, loop inside R is not adequate but single rules with expand comand would be sufficient in this pipeline.
- Let users to select if the analysis includes monomorphic sites, if appropriate.

1) angsd_sfs_mindepthind2.R
Description:
This scripts conduct 2) population-specific ANGSD with bamlist, 3) find intersect sites among different populations, 4) Re-run angsd on intersecting sites for each population, 5) estimate SFS using realSFS, 7) 2D SFS estimatino and Fst calculation using realSFS and associated functions, 8) clean up directories for log outputs.

Specification:
- I have commented out 2-4 since there is already a rules.angsd_sfs that runs on the intersecting sites found by angsd_intersect.smk. You can skip these sections, unless you find something you should add to the rules.angsd_sfs or the related bamlist rule.
- realSFS can be usable in the angsd environment (angsd.yaml). Please use the same environment.
- Please use the population definition in the config. Variablees like "continental" and "insular" are specific to my previous project.
- Let users to choose if SFS will be calculated using folded or unfolded spectrum. Set "fold" as default.
- Let users to choose if SFS will contain "linked" or "unlinked" sites. Set "linked" as default. If "unlinked" is set, then use the filtered sites found through rules.ngsld_global and its related rules.
- Let users to define maxiter, and tole. For thread, make it 10 under config/cluster.json. Set 16 for Memory size there.

2) diversity_indices_mindepthind2.R
Description:
This script reads outputs of realSFS and compare stats among different populations. It will save summary csv files, and also it conducts a statistical model to compare theta and tajima's D using lm, R package.

Specification:
- The script loads read_csv("/mnt/y/Studies/data/WGS/Laterallus_jamaicensis_coturniculus/contig_length.csv") that records contig sizes of each contig. However, this has not been calculated yet in this pipeline. Please obtain contig sizes of each contig and save it so that the pipeline will use this.
- Although it uses "group_2" to group samples, please stick with using the population defined in the config. Also, the script removes "kanto_w" here, please ignore this for the time being.
- Use R for the result organization and lm modelling. You can split the R scripts if needed.
- make a plot that compares Tajima's D and theta.
- Please output the summary of fitted models.
## Codex summary (appended by Codex on 2024-05-15)

### Files generated or updated
- `config/config_test.yaml`: added `sfs_analysis` section for site filters, folding options, monomorphic inclusion, and realSFS tuning.
- `config/cluster.json`: set dedicated resources for `angsd_sfs_group`, `realsfs_1d`, and `realsfs_2d` rules (10 threads / 16 GB as requested).
- `workflow/rules/common.smk`: normalised SFS configuration values, exposed shared flags (e.g., `SFS_NEEDS_UNLINKED_SITES`, `SFS_PAIR_LABELS`).
- `workflow/rules/angsd_sfs.smk`: replaced legacy single-rule setup with per-filter ANGSD, 1D/2D realSFS + Fst, and downstream summary rule.
- `workflow/rules/reference.smk`: added `export_contig_lengths` rule to emit contig lengths from the FASTA index.
- `workflow/Snakefile`: wired new SFS targets and contig-length export into `rule all`, conditionalised unlinked dependencies.
- `workflow/envs/plot.yaml`: included `r-broom` for tidy model outputs.
- `workflow/scripts/export_contig_lengths.py`: new helper script producing contig length CSV.
- `workflow/scripts/summarize_sfs.R`: new R script aggregating SFS spectra, θ/Tajima summaries, linear models, and plots.

### Description of changes
- Introduced configurable SFS analysis options (linked/unlinked, folded/unfolded, monomorphic inclusion) and passed them through shared globals.
- Refactored SFS workflow to emit SAF/MAF/BCF per population and to compute both 1D spectra/theta/Tajima and 2D spectra/Fst for all population pairs.
- Added pipeline step to export contig lengths from the reference `.fai` so downstream statistics can normalise by scaffold size.
- Implemented R-based summarisation that collates SFS outputs, fits `lm` models of θ and Tajima’s D vs population (with log-length covariate), and generates a comparison plot.

### Things to check before running & usage
- Ensure `sfs_analysis` settings in the active config reflect desired modes; set `site_filters: [linked, unlinked]` to trigger ngsLD-derived site usage.
- If `unlinked` spectra are requested, confirm the ngsLD branch has been executed so `results/unlinked_sites/...` is available (Snakemake handles ordering when full workflow is run).
- Run the SFS branch with e.g. `snakemake results/sfs/${output_prefix}/tajima_vs_theta.pdf --cores 16` (adjust cores ≥10 to match cluster profile for realSFS).
- Regenerate the full scenario via `snakemake --cores 30` to ensure new summaries integrate cleanly with other analyses.

### Additional notes
- `summarize_sfs.R` currently filters contigs <2 Mb; tweak `min_contig_length` in rule params if you need smaller scaffolds.
- The SFS summary plot facets by site filter and fold; expect empty facets when only one state is configured.
- Linear model output is written even when models cannot be fit (message recorded for transparency).

### Next steps
1. Review generated CSVs and `tajima_vs_theta.pdf` under `results/sfs/${output_prefix}` for sanity checks.
2. Update production configs to include desired fold/site combinations and re-run the workflow on real data.
3. Consider adding downstream reporting or archival steps once SFS outputs are validated.

_Changes appended at end of file by Codex._
