I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake. Please do this conversion only for the files specified below, since some of the scripts are already converted and integrated into this pipeline.

0) General
Descriptions & Specifications:
- These codes were developed in different projects. Please adapt them to the present snakemake pipelines on config parameters, tsv files to be read, directories, etc.
- cluster.script in the old script intends to convert shell commands to PBS cluster. Please ignore this step since snakemake is intended to be ran in PBS cluster in most cases (except when environment == "local").
- They are ran under environment == "cluster"
- You may add columns to any already-present tsv files to retrieve metadata of samples etc.
- Please adjust styles to snakemake. For example, loop inside R is not adequate but single rules with expand comand would be sufficient in this pipeline.

1) angsd_snapp_mindepthind2.R
This script prepares and runs SNAPP analysis based on ANGSD runs calling SNPs. The actual flow is as following.
- Select best 4 samples from each population, scoring the lowest missingness in the previous angsd global run
- Prepare a bamlist for the samples of all the populations (i.e., if four populations defined, 4*4 = 16 samples)
- Conduct ANGSD using the prepared bamlist on the intersecting sites defined in the angsd_intersect. 
- Obtain unlinked sites from the resultant ANGSD, based on rules from rules.ngsld_global ~ rules.angsd_filtered_sites_index (Not run the new ngsld)
- Remove excess heterozygosity sites using bcftools and convert vcf to phylip using vcf2phylip
- Create supplemental files in R and output an xml file using snapp_prep.rb.
- Run SNAPP analysis using an additional package of BEAST2 (outside the R script)
- Conduct treeannotator to get a summary tree
- Draw trees using ggtree

Specification
Please consider the following to implement this R script functionality to snakemake
- For ANGSD parameters, users can define them in the config like other ANGSD runs
- vcf2phylip and snapp_prep.rb are preapred under scripts
- SNAPP application is normally installed using "/path/to/beast/bin/packagemanager -add SNAPP". However, packagemanager cannot be used in the PBS cluster due to internet restrictions. Please consider any other ways to use SNAPP in the PBS cluster, or otherwise, let users run this locally.
- In the R script, runs under two conditions, "call/out1" and "call2/out2", were conducted. This is because they use different constraints. For the snakemake analysis, let users define one or multiple constraints in the config file (in the R script, the part stating lognormal(X,V,Y)), while type ("crown") and taxon (comma-delimited list of taxon) would be shared between different runs (but let users to define in the config file). Name different runs, like "run1", "run2" and make them run in parallel. 
- BEAST2 will be usable under conda env. Please create a yaml for this.
- For tree plot coloration, first forget about the color order of the populations. Label tips by outgroup names (like treemix) and population names (defined in the config).
- For the momment, ignore the part where treeannotator and ggtree were used to summarize trees and draw.  


### Codex Update – SNAPP Integration (appended by Codex)
1) **New/updated files**
   - Updated: `config/config_test.yaml`, `workflow/Snakefile`, `workflow/rules/common.smk`
   - New: `workflow/rules/snapp.smk`, `workflow/scripts/snapp_select_samples.py`, `workflow/scripts/snapp_extract_unlinked_sites.py`, `workflow/scripts/snapp_prepare_metadata.py`, `workflow/scripts/snapp_write_constraints.py`, `workflow/envs/snapp_prep.yaml`, `workflow/envs/beast_snapp.yaml`
2) **Change summary**
   - Added SNAPP-specific ANGSD parameters, sample selection knobs, and constraint definitions in config; exposed new runtime settings in `common.smk`
   - Wired a dedicated SNAPP rule chain (selection → ANGSD → filtering → phylip/XML prep → optional BEAST run) and registered targets with the main Snakefile
   - Implemented helper scripts for missingness ranking, unlinked extraction, metadata/constraint prep, plus minimal Conda envs for `snapp_prep.rb` and BEAST2/SNAPP
3) **Usage / things to check before running**
   - Update `config.snapp.enabled` to `true`, review per-population limits/aliases, and set constraint distributions you intend to run (e.g., `run1`, `run2`)
   - Ensure BEAST2 has the SNAPP plugin installed in the execution environment; enable `config.snapp.run_beast` only where BEAST can run (cluster queue vs. local)
   - To generate inputs without BEAST, run targets like `snakemake results/snapp/${output_prefix}/alignment/snapp.min${min_samples_locus}.phy`; include `--resources` as needed
4) **Additional notes**
   - Tree summarisation/visualisation (treeannotator, ggtree) remains to be ported
   - Outgroup inclusion defaults to sliced BAMs; adjust `snapp.include_outgroups` or `snapp.outgroup_samples` if a different panel is needed

---

### Additional configuration notes (2024-05-26)
- `snapp.outgroup_samples` must list the outgroup sample IDs (matching the `sample_id` column in `data/outgroups.tsv`) you want to add to the SNAPP bamlist. Leave it empty to omit outgroups; populate it with IDs such as `Ggallinago_01` when you need them included.
- `snapp.population_aliases` lets you remap ingroup population labels from the sample metadata (`group_col`) to the names SNAPP/BEAST should use. Provide a mapping like `megala_e: G_megala_E` when the analysis or constraints expect alternate labels.
- `snapp.outgroup_aliases` performs the same remapping for the outgroup entries defined in `snapp.outgroup_samples`, giving you control over the outgroup names that appear in metadata and exported files.
- Be sure to fill in both alias dictionaries (and `outgroup_samples`) before running the SNAPP workflow so downstream scripts can emit the desired identifiers.

### SNAPP XML-only workflow (2025-11-13)
- Repeated BEAST failures on the cluster were traced to BEAST 2.6.x lacking the SNAPP plugin and the cluster’s inability to run `packagemanager -add SNAPP` (no outbound network, Singularity user-dir not writable). Even after bundling a SIF, the runtime inside the cluster could not see `snap.Data`, so the job always died at XML parse time.
- To avoid that bottleneck, the pipeline now stops after XML preparation. Rule `snapp_prep_xml` writes `results/snapp/<profile>/runs/<run>/snapp.xml` plus per-run constraints. The `snapp.log_prefix` config key controls the relative `-prefix` passed to `snapp_prep` (defaults to `snapp`), so BEAST will emit `<prefix>.log`/`.trees` alongside the XML.
- There is no longer a `snapp_run_beast` rule, no `snapp.run_beast` toggle, and no Singularity/conda BEAST settings. Users must copy the XML (and constraints) to an environment where BEAST 2.7.x + SNAPP is installed (e.g., local macOS GUI), then run `beast -threads <n> snapp.xml`. TreeAnnotator/ggtree post-processing remains manual.
- When documenting SNAPP runs, note that the pipeline still automates sample selection, ANGSD calls, unlinked filtering, metadata prep, and XML generation; only the actual BEAST MCMC + summaries are out-of-band due to the plugin deployment limitations described above.
