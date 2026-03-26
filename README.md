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
git fetch upstream
git switch main
git pull origin main
git merge upstream/main
git push origin main
```
If you want to replace the upstream to other test branches, then
```bash
git fetch origin
git fetch upstream
git switch main
git branch --set-upstream-to=origin/main
git merge upstream/{test_branch_name}
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
The upstream repository keeps example/default files under:

- `templates/data/`
- `templates/config/`

These are reference templates, not active run files.

### Which files should you edit?

| Type | Use these files in your project | Notes |
| --- | --- | --- |
| Input metadata | `data/samples.tsv`, `data/outgroup.tsv`, `data/references.tsv` | Project-specific working files |
| Run config | `config/*.yaml` | Copy from `templates/config/` and edit |
| Launchers | `run_pipeline_local.sh`, `run_pipeline_cluster.sh` | Adjust to your environment as needed |
| Optional staging helper | `run_ssd_mirroring.sh` | Use only if you stage HDD data onto SSD locally |

This separation keeps upstream template updates distinct from project-specific metadata and configs.

> [!NOTE]
> Files under `config/validation_*.yaml` are developer validation fixtures for testing the pipeline itself. They are not intended to be copied as project configs.

## Step 5: Prepare your dataset
Refer to `templates/data/` when creating your project metadata under `data/`.

### Required project metadata

| File | Purpose | Key points |
| --- | --- | --- |
| `data/references.tsv` | Reference genome metadata | Use a closely related genome from NCBI Genome. `accession` is typically a `GCA...` accession. `name` should match the reference FASTA filename, usually ending in `genomic.fna.gz`. Add a `url` column only if you use `reference.download_method: wget`. |
| `data/samples.tsv` | Ingroup samples | `sample` is the downstream analysis ID. You may add additional population columns such as `pop2`, `pop3`, etc. |
| `data/outgroup.tsv` | Outgroup metadata / SRA accessions | If outgroups are not used, keep this file as a header-only TSV rather than a truly empty file. |

### Ingroup FASTQ naming and location
If all ingroup FASTQs live in one directory and follow a uniform naming scheme, you can still use `reads.ingroup_dir`.

If different samples live in different directories or use different naming rules, set per-sample columns in `samples.tsv` such as:

- `fastq_dir`
- `fastq_prefix`
- `fastq_r1_suffix`
- `fastq_r2_suffix`
- `fastq_extension`

The pipeline resolves ingroup FASTQ paths as:

```text
<fastq_dir>/<fastq_prefix><basename><fastq_r1_suffix><fastq_extension>
<fastq_dir>/<fastq_prefix><basename><fastq_r2_suffix><fastq_extension>
```

If preprocessing outputs should keep the original raw-data basename, add a separate column such as `basename` and point `reads.ingroup_metadata.basename_col` to it. Downstream analyses still use `sample`, while trimmed FASTQs and reusable BAMs follow `basename`.

### Adapters and reusable storage

1. Put adapter-related files under `data/adapters/`. See `templates/data/adapters/` for examples.
2. If you want reusable assets outside the repository, set `storage.*` in your config. These control where the pipeline stores:
   - references
   - downloaded / merged outgroup FASTQs
   - trimmed reads
   - QC outputs
   - temporary mapping files
   - filtered long reads
   - reusable BAMs

Downstream analysis outputs remain under the project-local `results/` and `figures/` directories.

## Step 5: Set up your own config files
Put project configs under `config/*.yaml`. Start from `templates/config/defaults_local.yaml` or `templates/config/defaults_cluster.yaml`, then edit the copy for your project.

### General guidance

| Recommendation | Why |
| --- | --- |
| Keep project configs under `config/` | Separates project state from upstream templates |
| Create separate configs for separate runs | Makes parameter changes and reruns traceable |
| Set a distinct `project.output_prefix` for each run | Keeps outputs separated under `results/` and `figures/` |

### Storage roots and shared layouts
The `storage` section supports two common patterns.

| Pattern | When to use it | Typical keys |
| --- | --- | --- |
| Explicit absolute directories | You want fine control per asset type | `storage.reference.dir`, `storage.bam.ingroup_dir`, `storage.bam.outgroup_dir`, etc. |
| Shared relative layout under different roots | Local and cluster share the same internal directory structure but different root prefixes | `storage.roots.local`, `storage.roots.cluster`, `storage.shared.*` |

Example:

```yaml
storage:
  roots:
    local: /Volumes/OWCEnvoyProFX/mirror_data
    cluster: /lfs/aokid/data
  shared:
    reference:
      dir: source/reference/GCA_020086605.1
    bam:
      ingroup_dir: ngs/derived/project_a/v1/bam
      outgroup_dir: ngs/derived/project_a/v1/bam/outgroups
```

With this layout, the same relative subpath is resolved beneath different roots depending on the active environment.

### Notes on `storage.shared.*`

- `storage.shared.*` entries are treated as relative subdirectory layouts.
- If `storage.roots` is null, `storage.shared.*` is ignored.
- In that case, the workflow falls back to explicit `storage.*` keys or built-in defaults.

### Ingroup FASTQ metadata
The `reads.ingroup_metadata` section controls how per-sample FASTQ paths are resolved from `samples.tsv`.

| Config key | Default TSV column |
| --- | --- |
| `basename_col` | `sample` |
| `dir_col` | `fastq_dir` |
| `prefix_col` | `fastq_prefix` |
| `r1_suffix_col` | `fastq_r1_suffix` |
| `r2_suffix_col` | `fastq_r2_suffix` |
| `extension_col` | `fastq_extension` |

If those per-sample columns are absent, the workflow falls back to `reads.ingroup_dir` plus `_1`, `_2`, and `.fastq.gz`.

### Transfer mirror root
`transfer.cluster_storage_root` controls the root directory used for mirrored external paths on the cluster.

Example:

```yaml
paths:
  transfer:
    cluster_storage_root: /lfs/aokid/data
```

With:

```yaml
storage:
  roots:
    local: /Volumes/SSD
    cluster: /lfs/aokid/data
```

a local path such as `/Volumes/SSD/project_a/reference` is mirrored on the cluster as `/lfs/aokid/data/project_a/reference`.

### Analysis toggles
Each analysis block has an `enabled: true|false` switch. Use these to turn analyses on or off without editing the Snakefile.

When no outgroups are used:

- set `include_outgroups: false` for analyses such as global ANGSD, RAxML, TreeMix, ngsDist, and SNAPP
- disable ABBABABA2 unless an outgroup panel is provided

## Step 6: Adjust run_pipeline_cluster.sh/run_pipeline_local.sh
These two bash scripts are prepared to run Snakemake entirely either in a cluster or a local environment. While they make it easy to run from the beginning to the end of this pipeline, you may want to change some details sometimes. Please adjust the contents, especially options of `snakemake` to your own purposes.
For cluster execution, the current Slurm setup uses `profile/default/config.yaml` together with rule-level `threads` and `resources`. `templates/config/cluster.json` is kept only as a legacy PBS reference and is not used by `run_pipeline_cluster.sh`.

If you keep cold storage on HDD but want better local read/write performance on SSD, you may stage those files with `run_ssd_mirroring.sh` before running the local workflow. The script now reads a TSV manifest, defaults to dry-run preview mode, performs the actual copy only when called with `--apply`, and writes timestamped rsync logs next to the manifest by default.

Recommended locations:

| File | Purpose |
| --- | --- |
| `data/ssd_mirroring.tsv` | Project-specific working manifest that you edit |
| `templates/data/ssd_mirroring.tsv` | Example template showing the expected columns |

The TSV should contain:

| Column | Meaning |
| --- | --- |
| `cold_path` | Source path on cold storage |
| `ssd_path` | Destination path on SSD |

Example:
```bash
bash run_ssd_mirroring.sh data/ssd_mirroring.tsv
bash run_ssd_mirroring.sh --apply data/ssd_mirroring.tsv
bash run_transfer_verification.sh data/ssd_mirroring.tsv
```

When using this staging pattern, point your local config paths such as `storage.roots.local` or `reads.ingroup_dir` to the SSD-side paths that the pipeline should actually read and write.

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
The workflow generates transfer helpers under `transfer/<output_prefix>/`.

### Generated files

| File | Purpose |
| --- | --- |
| `external_paths.tsv` | Manifest of external files/directories and their mirrored cluster paths |
| `sync_local_external_to_cluster.sh` | Upload local external files to the cluster mirror root |
| `sync_cluster_external_to_local.sh` | Copy mirrored cluster files back to the original local external paths |
| `project_output_paths.tsv` | Manifest of project-level output directories mirrored from the cluster repo |
| `sync_cluster_project_outputs_to_local.sh` | Pull `results/`, `figures/`, and `logs/` from the cluster project into the local project |
| `cold_storage_paths.tsv` | Manifest of external SSD paths that can be back-mirrored into cold storage |

Typical use is as follows.
```bash
# SSD -> Cluster (raw shared files)
bash transfer/defaults/sync_local_external_to_cluster.sh username@cluster.example.org
# Cluster -> SSD (derived shared files)
bash transfer/defaults/sync_cluster_external_to_local.sh username@cluster.example.org
# Cluster -> SSD (results, figures, logs)
bash transfer/defaults/sync_cluster_project_outputs_to_local.sh username@cluster.example.org /lfs/aokid/project_a
# SSD -> HDD
bash run_cold_storage_mirroring.sh transfer/defaults/cold_storage_paths.tsv /Volumes/cold_storage/project_a
# Verify SSD -> HDD state
bash run_transfer_verification.sh transfer/defaults/cold_storage_paths.tsv /Volumes/cold_storage/project_a
```

Recommended config:

```yaml
paths:
  transfer:
    cluster_storage_root: /lfs/aokid/data
    cluster_project_dir: /lfs/aokid/project_a
```

### Which machine should run each script?

| Script | Where to run it | Why |
| --- | --- | --- |
| `sync_local_external_to_cluster.sh` | Local machine | The script reads local files and pushes them to the cluster |
| `sync_cluster_external_to_local.sh` | Local machine | The script pulls from the cluster back to the original local storage paths |
| `sync_cluster_project_outputs_to_local.sh` | Local machine | The script pulls `results/`, `figures/`, and `logs/` from the cluster project into the local project |
| `run_cold_storage_mirroring.sh` | Local machine | The standalone utility reads `cold_storage_paths.tsv` and back-mirrors SSD-resident external files into cold storage |
| `run_transfer_verification.sh` | Local machine | The standalone verifier compares manifest source and destination paths and writes a TSV status report |

> [!IMPORTANT]
> `sync_cluster_external_to_local.sh` should normally be generated and run from the local analysis if you want the embedded destination paths to point to your SSD, such as `/Volumes/OWCEnvoyProFX/...`. If you generate it on the cluster, its baked-in “local” paths will instead reflect cluster-resolved paths such as `/lfs/aokid/...`.

### Optional environment variables

| Variable | Meaning |
| --- | --- |
| `CLUSTER_PORT` | SSH port for `rsync` / `ssh` |
| `CLUSTER_STORAGE_ROOT` | Override the mirror root without editing config |
| `CLUSTER_PROJECT_ROOT` | Override the cluster-side project directory for project output download |

# Tips
## Dag and dryrun with "mock" option
```bash
snakemake -n -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true
snakemake --dag -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true | dot -Tpdf > dag.pdf
```
