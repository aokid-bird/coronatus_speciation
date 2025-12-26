I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake. Please do this conversion only for the files specified below, since some of the scripts are already converted and integrated into this pipeline.

0) General
Descriptions & Specifications:
- These codes were developed in different projects. Please adapt them to the present snakemake pipelines on config parameters, tsv files to be read, directories, population/sample names etc.
- cluster.script in the old script intends to convert shell commands to PBS cluster. Please ignore this step since snakemake is intended to be ran in PBS cluster in most cases (except when environment == "local").
- They are ran under environment == "cluster"
- You may add columns to any already-present tsv files to retrieve metadata of samples etc.
- Please adjust styles to snakemake. For example, loop inside R is not adequate but single rules with expand comand would be sufficient in this pipeline.
- Let users to select if the analysis includes monomorphic sites, if appropriate.
- Parameters inside the R script, especially population names and individual names, are specific to the previous projects. The pipeline should be generalized to use many purposes. These values should be read from the configs and tsv files. If any specific parts are detected, please learn things from the present pipeline to understand what may indicate the populations/individuals. If anything is unclear for you, please ask me before changes are made.

1) angsd_abbababa.R
The script generates a bamlist file, calculates parameters needed for abbababa2 analysis (sizefile, sample.name) and save them as files, and run abbababa2. Please do the following to implement the functionality of this script to the snakemake pipeline.
- To generate the bamlist, let the users to choose which outgroups to include (sample_id).
- Please use the intersecting sites inferred from rules around angsd_intersect.smk to do abbababa and refer to it in -sites and -rf
- let users to define the label of outgroup. If NULL, then convert "taxon" from the outgroups.tsv to any appropriate label(s).
- Use the codes under "3) D-stat" of the R script to calculate D-stats. code/estAvgError.R in the script refers to workflow/cripts/estAvgError.R. Please do not make any change inside estAvgError.R. You may separate the part where estAvgError.R is used and where dataframe organization is done.

2) figures_abbababa_mindepthind2.R
- Use khroma for the palette. Please refer to any other snakemake rules to implement khroma because this cannot be used under a conda env.
- Use population labels from the config file.

