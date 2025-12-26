## Step 0: Prepare your environments
Install miniconda
https://www.anaconda.com/docs/getting-started/miniconda/install#macos

## Step 1: Git clone the repository
```bash
# store the token to credential helper
git config --global credential.helper store
# for the first time
git clone https://github.com/aokid-bird/gbs_pipeline.git
# to pull the latest version (main)
git pull origin main
```

## Step 2: Set up Snakemake environment
In the current snakemake.yaml, the snakemake environment (bioinfo_pipeline) will be created with python 3.10 and Snakemake 7.32.4. This will allow you to use f-string in the Snakemake rules.
```bash
conda env create -f snakemake.yaml
conda activate bioinfo_pipeline
```

## Step 3: Prepare your dataset
Please refer to any default files or `templates/data` to prepare your own input metadata.
1. Find a suitable reference genome from NCBI.Genome. Prepare your reference list `data/reference.tsv`. Look up Genome of NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/) and use closely related species. `accession` column is something with "GCA", and `name` should include something ending with "genomic.fna.gz" which you can usually find in the ftp tab of the genome. If you intend to use `config.reference_download_method = wget` option, then create a column called `url` to retrieve the data. Delete the `url` column if you use the `dataset` method. This file should be referred to at `config.references_tsv`.
2. Prepare your sample lists `data/sample.tsv`. If there are more than one population definition, then create new columns like `pop2`, `pop3`..., which will be referred to in the `config/config.yaml`. Note that the sample name will be used to identify the fastq files, bam files, etc, so unify the names between fastq files and this lists. This file should be referred to at the config.samples.
3. Prepare your outgroup list `data/outgroup.tsv`. Look up SRA of NCBI (https://www.ncbi.nlm.nih.gov/sra) with your keywords, like genus name. Fill each column. For example, normally, `sample_id` starts with "SAMN", `srr_id` starts with "SRR".
4. Put your data in `data/raw`. Make sure that the names of fastq files match with the name lists in the sample.tsv.
5. Put your adapter sequence fasta in `data/adapters`. Please refer to `templates/data/adapters` for an example.

## Step 4: Set up your own config files
The config file is the most important part of the analysis where you define specific parameters for each analysis. Put your config files in config/XXX.yaml. Please refer to templates/config for an example. The `config/defaults_cluster.yaml` or `config/defaults_local.yaml` will guide you to make your own. 
- Create your own config file for each specific tasks and runs. Please make separate config files when you want to change parameters, change population definitions, etc. Then, set different names to `config.output_prefix` whose unique directory will be created under `results/ANALYSISNAME/OUTPUT_PREFIX`.
- The pipeline is still under development, and therefore, flexibility is still low. In the future, you may choose which analyses you want to do by turning "enabled: true" on. In the present version, however, you may need to go through all or most of the analyses.

## Step 5: Adjust run_pipeline_cluster.sh/run_pipeline_local.sh
These two bash scripts are prepared to run Snakemake entirely either in a cluster or a local environment. While they make it easy to run from the beginning to the end of this pipeline, you may want to change some details sometimes. Please adjust the contents, especially options of `snakemake` to your own purposes.

## Step 6: Analysis environment
Once creating your snakemake environment `bioinfo_pipeline`, other analyses environment will be created inside the pipeline (for cases of conda) or before the pipeline initiates (for cases of singularity using the run_pipeline.sh). Therefore, you do not have to create your own analysis environment. 

## Step 7: Place software bin files and FASTQ, Reference, SRA-FASTQ and other huge files to the cluseter by using FTP transfer.
Because some files are too huge to place via github, you need to manually transfer files like .bin/.sif files (mainly softwares) and raw data (fastq, .fasta and their index files) to your cluster on your own. You normally use FTP and any other protocols to transfer between your local and cluster.

## Step X: Run your snakemake
Once your bash scripts are prepared, inside the pipeline root directory, simply type the following.
```bash
# inside your local environment (to download genome/SRA)
bash run_pipeline_local.sh
# inside your cluster environment (for the main run)
bash run_pipeline_cluster.sh
```

# Tips
## Dag and dryrun with "mock" option
```bash
snakemake -n -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true
snakemake --dag -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true | dot -Tpdf > dag.pdf
```
