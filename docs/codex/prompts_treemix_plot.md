I am developing a snakemake pipeline for genomics biofinformatics and analyses. Meanwhile, in the folder "codes_from_old_projects", I included some R scripts that I used for similar purposes work similarly to the snakemake but codes are lengthy and redundant. Things included in this folder are some of the scripts I used. Can you convert them either to snakemake rules, config parameters, parameter manipulations in common.smk, to be integrated to this snakemake pipeline? I am providing brief descriptions to general and each R scripts and some specifications you may follow to adapt them to snakemake. Please do this conversion only for the files specified below, since some of the scripts are already converted and integrated into this pipeline.

0) General
Descriptions & Specifications:
- These codes were developed in different projects. Please adapt them to the present snakemake pipelines on config parameters, tsv files to be read, directories, etc.
- cluster.script in the old script intends to convert shell commands to PBS cluster. Please ignore this step since snakemake is intended to be ran in PBS cluster in most cases (except when environment == "local").
- They are ran under environment == "cluster"
- You may add columns to any already-present tsv files to retrieve metadata of samples etc.
- Please adjust styles to snakemake. For example, loop inside R is not adequate but single rules with expand comand would be sufficient in this pipeline.


1) codes_from_old_projects/figures_2_rlrmv_mindepthind2_treemix.R (= previous file in the following)
- This R script includes codes to plot treemix or evaluation of treemix results. 
- The optM packages, scripts/treemix_plotting_funcs.R (equivalent with "~/treemix-1.13/src/plotting_funcs.R"), scripts/ggplot_treemix.R will be used inside by sourcing scripts or reading R packages. 
- The script is written for the previous projects, and some of the variables are needed to be generalized. 
a. Outgroups should be defined in the config, but the previous file is set to be "Laterallus jamaicensis" or "Laja" as its abbreviation.
b. bamlist.df in the previous file should now be read from the one of the bamlist rules and the tsv file set in the config.
c. Number of replicates (in the previous file, 10) and the number of migration edges tested (in the previous file, 6) should be generalized by receiving variables from snakemake
d. Please use khroma as a color R package. But in this pipeline, use "bright" palette in the original order, rather than specifying the manual order.
- Although there are redundant codes, please let them as is. Simply embed the plotting codes inside the snakemake pipeline without your original modifications. If need some, then please ask me before implementing. 
- For saving plot, adjust the codes to adapt them to the cluster environment. This is because ggsave may not work in the cluster environment. 
