## Step 0: Prepare your environments
Install miniconda
https://www.anaconda.com/docs/getting-started/miniconda/install#macos


## Step 1: Create your project directory in the local
```bash
mkdir -p ~/proj_name
```
Open this project directory through your VScode.

## Step 2: Git clone the repository
```bash
# enter the project dir
cd ~/peoj_name
# store the token to credential helper
git config --global credential.helper store
# for the first time
git clone https://github.com/aokid-bird/gbs_pipeline.git .
# to pull the latest version (main)
git pull origin main
```
Because this pipeline will be used in many situations and these runs will require furhter updates, this repository should be used as "submodule" of each project.

> [!CAUTION]
> project *gallinag_phylogeography* and *coronatus_speciation* are not updated as submodules yet (2025/12/26). When rerunning the analyses in these projects, please update your directories by making new directories and downloading this pipeline as submodules, and transfer some input/config files already generated in the older project directories. 

## Step 3: Set the template as Upstream and create the projectory-specific repos
```bash
# rename the template
git remote rename origin upstream
# create your github repo on the web browser

# push the project repo to github
git remote add origin git@github.com:YOUR_GITHUB_USERNAME/proj_name.git
# then push the repo from the local to remote
# -u option is to declare that upstream is added 
git push -u origin main
```
This will create a repo for your own project, while this is based on the template gbs_pipeline.
The following results of `git remote -v` is correct to push updates to your project repo (not the template)
```bash
git remote -v
> origin  git@github.com:YOUR_GITHUB_USERNAME/proj_name.git (fetch)
> origin  git@github.com:YOUR_GITHUB_USERNAME/proj_name.git (push)
> upstream        https://github.com/aokid-bird/gbs_pipeline.git (fetch)
> upstream        https://github.com/aokid-bird/gbs_pipeline.git (push)
```

## Step 3.5: Update your template
When the template gbs_pipeline is updated by the creater, you may update your project by updating the template.
```bash
git checkout main
git fetch upstream
git merge upstream/main
git push origin main
```

## Step 4: Set up Snakemake environment
In the current snakemake.yaml, the snakemake environment (bioinfo_pipeline) will be created with python 3.10 and Snakemake 7.32.4. This will allow you to use f-string in the Snakemake rules.
```bash
conda env create -f snakemake.yaml
conda activate bioinfo_pipeline
# install graphviz for dot used for drawing DAG
brew install graphviz
```

## Project files versus templates
The upstream pipeline repository keeps example/default files under:
- `templates/data/`
- `templates/config/`

These are reference templates, not the active files for a project run.

For each project repository, create and edit your working files under:
- `data/samples.tsv`
- `data/outgroup.tsv`
- `data/references.tsv`
- `config/*.yaml`
- `run_pipeline_local.sh`
- `run_pipeline_cluster.sh`

This keeps upstream template updates separate from project-specific metadata and config changes.

The files under `config/validation_*.yaml` in this pipeline repository are developer validation fixtures used for dry-run testing of the workflow itself. They are not intended to be copied as project configs.

## Step 5: Prepare your dataset
Please refer to `templates/data/` when creating your own project input metadata under `data/`.
1. Find a suitable reference genome from NCBI.Genome. Prepare your reference list `data/references.tsv`. Look up Genome of NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/) and use closely related species. `accession` column is something with "GCA", and `name` should include something ending with "genomic.fna.gz" which you can usually find in the ftp tab of the genome. If you intend to use `config.reference_download_method = wget` option, then create a column called `url` to retrieve the data. Delete the `url` column if you use the `dataset` method. This file should be referred to at `config.references_tsv`.
2. Prepare your sample list `data/samples.tsv`. If there are more than one population definition, then create new columns like `pop2`, `pop3`..., which will be referred to in the `config/config.yaml`. The `sample` column remains the canonical sample ID used for BAMs and downstream analyses.
For ingroup FASTQ discovery, you may now also define per-sample columns such as `fastq_dir`, `fastq_prefix`, `fastq_r1_suffix`, `fastq_r2_suffix`, and `fastq_extension`. The pipeline resolves each sample path as:
`<fastq_dir>/<fastq_prefix><sample><fastq_r1_suffix><fastq_extension>`
and
`<fastq_dir>/<fastq_prefix><sample><fastq_r2_suffix><fastq_extension>`.
This allows mixed storage locations and file naming rules within one run.
3. Prepare your outgroup list `data/outgroup.tsv`. Look up SRA of NCBI (https://www.ncbi.nlm.nih.gov/sra) with your keywords, like genus name. Fill each column. For example, normally, `sample_id` starts with "SAMN", `srr_id` starts with "SRR".
4. Put your data where you want. If all ingroup FASTQs follow a single shared directory and default naming rule, you may still use `reads.ingroup_dir`. If samples are distributed across multiple directories or use different naming conventions, define those per sample in `samples.tsv` and the pipeline will use those metadata instead.
5. Put your adapter sequence fasta in `data/adapters`. Please refer to `templates/data/adapters` for an example.
6. If you want reusable assets outside the repository, set the explicit `storage.*` keys in your config. These control where the pipeline stores reference files, downloaded outgroup FASTQs, merged FASTQs, trimmed reads, QC outputs, temporary mapping files, filtered long reads, and reusable BAMs. Downstream analysis outputs such as ANGSD, RAxML, TreeMix, SNAPP, and plots remain under the project-local `results/` and `figures/` directories.

## Step 5: Set up your own config files
The config file is the most important part of the analysis where you define specific parameters for each analysis. Put your project config files in `config/XXX.yaml`. Please refer to `templates/config/` for examples. The `templates/config/defaults_cluster.yaml` and `templates/config/defaults_local.yaml` files are reference templates and should be copied into project-specific config files before editing.
- Once you make your own config, put them under `config/`. 
- Create your own config file for each specific tasks and runs. Please make separate config files when you want to change parameters, change population definitions, etc. Then, set different names to `config.output_prefix` whose unique directory will be created under `results/ANALYSISNAME/OUTPUT_PREFIX`.
- The new `storage` section uses explicit keys. Leave them at their defaults to keep files inside the repository, or point them to absolute directories like `/Data/WGS/...` for reusable storage across analyses. If `storage.bam.ingroup_dir` or `storage.bam.outgroup_dir` is set, those BAMs are treated as reusable preprocessing outputs rather than scenario-specific files.
- If local and cluster use different root prefixes, you may instead set `storage.roots.local` and `storage.roots.cluster` together with the relative layout under `storage.shared.*`. This lets local paths such as `/Volumes/...` and cluster paths such as `/lfs/aokid/...` share the same internal directory structure without repeating every absolute path. Explicit per-directory keys like `storage.reference.dir` still override that root-based layout when needed.
- `storage.shared.*` entries are treated as relative subdirectory structures under the active root. If `storage.roots` is left null, those shared entries are ignored and the workflow instead uses explicit `storage.*` / `reads.ingroup_dir` values when provided, otherwise it falls back to the built-in in-repo defaults.
- The `reads.ingroup_metadata` section defines which `samples.tsv` columns should be used for ingroup FASTQ directory and naming resolution. The defaults expect columns named `fastq_dir`, `fastq_prefix`, `fastq_r1_suffix`, `fastq_r2_suffix`, and `fastq_extension`, with fallback to `reads.ingroup_dir` plus `_1`/`_2` and `.fastq.gz` when those per-sample columns are absent.
- The `transfer.cluster_storage_root` setting defines where external absolute paths are mirrored on the cluster. For example, a local external path `/Data/WGS/project_a/reference` is transferred to `/lfs/aokid/Data/WGS/project_a/reference` when `transfer.cluster_storage_root` is `/lfs/aokid`.
- Each analysis block now has an `enabled: true|false` switch in the config templates, so you can turn major downstream analyses on or off per run without editing the Snakefile.

## Step 6: Adjust run_pipeline_cluster.sh/run_pipeline_local.sh
These two bash scripts are prepared to run Snakemake entirely either in a cluster or a local environment. While they make it easy to run from the beginning to the end of this pipeline, you may want to change some details sometimes. Please adjust the contents, especially options of `snakemake` to your own purposes.
For cluster execution, the current Slurm setup uses `profile/default/config.yaml` together with rule-level `threads` and `resources`. `templates/config/cluster.json` is kept only as a legacy PBS reference and is not used by `run_pipeline_cluster.sh`.

## Step 7: Analysis environment
Once creating your snakemake environment `bioinfo_pipeline`, other analyses environment will be created inside the pipeline (for cases of conda) or before the pipeline initiates (for cases of singularity using the run_pipeline.sh). Therefore, you do not have to create your own analysis environment. 
The workflow now standardizes generic helper environments into:
- `workflow/envs/python_utils.yaml` for pandas-based Python helper scripts
- `workflow/envs/r_plotting.yaml` for generic R plotting and reporting scripts
Specialized environments such as TreeMix evaluation, RAxML tree plotting, ABBABABA2 error estimation, and SNAPP preparation remain separate.

## Step 8: Place software bin files and FASTQ, Reference, SRA-FASTQ and other huge files to the cluseter by using FTP transfer.
Because some files are too huge to place via github, you need to manually transfer files like .bin/.sif files (mainly softwares) and raw data (fastq, .fasta and their index files) to your cluster on your own. You normally use FTP and any other protocols to transfer between your local and cluster.
Each externalized storage directory also receives a `README.md` and `provenance.yaml` written by the workflow. These files record which pipeline wrote the directory, the project path, the `output_prefix` active at the time, and the type of reusable assets stored there.

The bin files should be treated as following, after relocating them to the cluster. Do this once in the project repository, and it will be fine.

```bash
chmod +x workflow/bin/glactools
```

## Step X: Run your snakemake
Once your bash scripts are prepared, inside the pipeline root directory, simply type the following.
```bash
# inside your local environment (to download genome/SRA)
bash run_pipeline_local.sh
# inside your cluster environment (for the main run)
bash run_pipeline_cluster.sh
```

To test the Slurm path directly without the wrapper:
```bash
snakemake --profile profile/default --snakefile workflow/Snakefile --configfile config/config.yaml --dryrun
```

## Transfer helper scripts
The workflow generates transfer helpers under `results/transfer/<output_prefix>/`:
- `external_paths.tsv`: manifest of external files/directories and their mirrored cluster paths
- `sync_local_external_to_cluster.sh`: upload local external files to the cluster mirror root
- `sync_cluster_external_to_local.sh`: copy mirrored cluster files back to the original local external paths

Typical use:
```bash
bash results/transfer/defaults/sync_local_external_to_cluster.sh username@cluster.example.org
bash results/transfer/defaults/sync_cluster_external_to_local.sh username@cluster.example.org
```

Optional environment variables:
- `CLUSTER_PORT`: SSH port for `rsync`/`ssh`
- `CLUSTER_STORAGE_ROOT`: override the mirror root without editing config

# Tips
## Dag and dryrun with "mock" option
```bash
snakemake -n -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true
snakemake --dag -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true | dot -Tpdf > dag.pdf
```
