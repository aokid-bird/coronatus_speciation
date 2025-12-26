## Step 0: Prepare your environments
qInstall miniconda
https://www.anaconda.com/docs/getting-started/miniconda/install#macos

## Step 1: Git clone the repository
```bash
# store the token to credential helper
git config --global credential.helper store
# for the first time
git clone https://github.com/aokid-bird/gallinago_phylogeography.git
# to pull the latest version (main)
git pull origin main
# to pull the test version
git pull origin test
```

## Step 2: Set up Snakemake environment
In the current snakemake.yaml, the snakemake environment (bioinfo_pipeline) will be created with python 3.10 and Snakemake 7.32.4. This will allow you to use f-string in the Snakemake rules.
```bash
conda env create -f snakemake.yaml
conda activate bioinfo_pipeline
```

## Step 3: Set up conda analysis environment
Step 2 created a conda environment for Snakemake to run. It will call many different analyses that will require specific environments, which will be created by either conda (this step) or [Step 4](#step-4-obtain-sif-image-from-docker). Conda environments are defined under ```workflow/envs/``` and you can create an environment based on the yaml files. For most cases, these envrionments are created automatically while running the pipelines.
```bash
conda env create -f workflow/envs/filename_of_env.yaml
```

## Step 4: Obtain .sif image from docker
Some analyses require manual instalation of tools, which is difficult in the supercomputer. Furthermore, it is difficult to prepare exactly the same environment by the manual installation. Therefore, we use singularity. For most cases, singularity pre-made images are downloaded and saved in /home/user/envs/ in the following steps. This is also done in the bash script "run_pipeline.sh", so it does not normally require manual downloads.
```bash
mkdir home/{yourusername}/envs # if not created yet
cd home/{yourusername}/envs
singularity build filename.sif docker://address/to/image
``` 
The bash scripts are prepared to generate the singularity .sif images, being stored in ```workflow/containers/{analysis_name}/build.sh``` and their metadata are stored in ```workflow/containers/containers.yaml```. 
Simply execute
```bash
bash workflow/containers/{analysis_name}/build.sh
```

## Step 5: Create a config file
The config file is the most important part of the analysis where you define specific parameters for each analysis. The pipeline is still under development, and therefore, flexibility is still low. In the future, you may choose which analyses you want to do by turning "enabled: true" on. In the present version, however, you may need to go through all or most of the analyses.

1. Set the reference genome information. If config.reference.fasta is null, then it will read data/reference.tsv to retrieve reference fasta. Note that this functionality is available only in the local environment (config.environment = local). Find a suitable reference genome from NCBI.Genome. Prepare your reference list "reference.tsv". Look up Genome of NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/). Use closely related species name to look up. "accession" column is something with "GCA", and "name" should include something ending with "genomic.fna.gz" which you can usually find in the ftp tab of the genome. If you intend to use "config.reference_download_method = wget" option, then fill url with the url to retrieve the data. This file should be referred to at config.references_tsv
2. Prepare your sample lists "sample.tsv". If there are more than one population definition, then create new columns like "pop2", "pop3"..., which will be referred to in the config.yaml. Note that the sample name will be used to identify the fastq files, bam files, etc. This file should be referred to at the config.samples.
3. Prepare your outgroup list "outgroup.tsv". Look up SRA of NCBI (https://www.ncbi.nlm.nih.gov/sra) with your keywords, like genus name. Fill each column. For example, normally, "sample_id" starts with "SAMN", "srr_id" starts with "SRR".
4. Create your own config file for each specific tasks and runs. Please make separate config files when you want to change parameters, change population definitions, etc. Then, set different names to "output_prefix" whose unique directory will be created under results/ANALYSISNAME/OUTPUT_PREFIX.

## Dag and dryrun with "mock" option
```bash
snakemake -n -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true
snakemake --dag -s workflow/Snakefile -j 1 -p --config dryrun_mock_reference=true | dot -Tpdf > dag.pdf
```
