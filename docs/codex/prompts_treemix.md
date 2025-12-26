I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake. Please do this conversion only for the files specified below, since some of the scripts are already converted and integrated into this pipeline.

0) General
Descriptions & Specifications:
- These codes were developed in different projects. Please adapt them to the present snakemake pipelines on config parameters, tsv files to be read, directories, etc.
- cluster.script in the old script intends to convert shell commands to PBS cluster. Please ignore this step since snakemake is intended to be ran in PBS cluster in most cases (except when environment == "local").
- They are ran under environment == "cluster"
- You may add columns to any already-present tsv files to retrieve metadata of samples etc.
- Please adjust styles to snakemake. For example, loop inside R is not adequate but single rules with expand comand would be sufficient in this pipeline.

1) angsd_treemix_mindepthind2_minind0.5.R
Descriptions
This script conducts 2) ANGSD and output a bcf file, 3) convert bcf to acf, while specifying populations, then generate a treemix input using glactools, 3-5) calculate block size parameter used in treemix, using my original script called lddecay_blocksizefinder.R, 4) Run treemix, 5) evaluate models using plotting_func.R implemented in treemix, R optM package, and my own functions
Specifications
- Binary version of glactools will be used. This binary file will be transferred to the PBS cluster manually and stored in workflow/bin. Since .gitignore includes workflow/bin, this bin file will not be pushed to git, but makesure the use of glactools from this directory. 
- This script runs multiple sets of treemix for different data subsets. Please ignore this in the snakemake pipeline because different subsets should be ran using the same pipeline with different config settings
- I have prepared a smk file called "angsd_treemix.smk". Please modify this for pipeline adaptation.
- I want to include ougroups to this analysis. This R script does this, while outgroup implementation to the current version of the snake pipeline has not been done yet. Essentially, outgroups are defined in the outgroups.tsv.
- Let users to choose whether to include outgroups and whether to conduct per-population treemix or per-individual treemix, like this R script does in the loop.
- Please use lddecay_blocksizefinder.R or its modified version. I rather stick to use R for this analysis.
- "~/treemix-1.13/src/plotting_funcs.R" is provided by treemix. Should the pipeline include this?

Please make a document called "docs/codex/treemix_angsd.md" that include 1) my prompts, 2) your responses and descriptions on changes, and everything you have done in this chat. Then please update docs/codex/index.md.
