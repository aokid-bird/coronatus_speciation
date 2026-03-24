I am planning to refactor the present pipeline. This pipeline was created by both the users and codex that worked to convert my previous R pipelines to adapt them to snakemake. This was done step-by-step through different rules and methods. These backgrounds resulted in the fluctuation in the coding styles. In this refactoring, please unify the coding styles across the pipeline. In order to accomplish this, please follow the following points.

Before doing anything, please make agenda to meet the following. If necessarry, you can unite several prompts into one for better flows, or you can also prompts to different agenda if you think it appropriate. If anything is stated unclearly, ask the user to clarify it.

1. Please group the rules into functional means. This has been basically done by making different rules/.smk files. However, it is possible that several functionally seggreated rules are lumped into one .smk file. Therefore, please evaluate the present grouping, and list groups with some explanations. Do this and let users to decide whether to accept, before advancing the refactoring any further.
2. Organize `templates/config/defaults_cluster.yaml` and `templates/config/defaults_local.yaml` to follow the grouping determined in 1. However, please let the following sections grouped regardless of the grouping determined in 1. Please stick with the levels of the config (like, config.level1.level2.level3...) if adequate, but if grouping is much better, you may change the levels and nestedness.
- Basic settings and output_prefix for the analysis set
- Directories for reference and bams and other general directories
- reference and sample information and tsv file paths
- population and outgroup information in general (exclude those specific to different anlayses)
- angsd arguments should nest different arguments for different angsd analyses. Please move angsd_args currently under the specific analysis configs to this too (e.g., abbababa2.angs_args -> angs_args.abbababa2)
- minIndRatio for angsd
3. Pleae add config.[somerules].enabled for each analysis grouped in 1. Accordingly, add the functionality of this to each rule to allow users to decide which analyses they will conduct.
4. Please add include_outgroups to config.ngsdist and modify the rules accordingly.
5. Many rules include the functionality to remove some samples/outgroups and/or set specific outgroup samples as "outgroups". However, the ways to accomplish this may differ among analyses. Please unify the methods. It is okay to have the methods for each analysis group, or if adequate and if making it readable codees, you may also physically unify them by making common scripts.
6. rule.all is messy. Organize them. 
7. There are many variables defined in the common.smk. However, many snakemake rules and Snakefile define variables too. This is complicated and difficult to follow where the variables are defined. Also, many helper functions scatter inside rules (e.g. angsd_sfs.smk, and others too). Ideally, snakemake may prefer simpler styles, for readability and trackability. Please standardize whether, how, when and what variables and helpers are created. If not necessary, then use simpler ways; e.g. use f-string just inside the rules.input, rules.output etc without using helpers/pre-defined variables. If, for flexibility, variables and helpers are essential, then use them but in readable and organized manners. 
8. Please make Snakefile simpler. 
9. Please add precise explanations to each rule at the beginning of the rule """/"""
10. Now figures are saved either to results or figures. Please standardize across rules to save figures under figures. Please follow raxml rules for the directory locations. Standardize figures to pdfs, not pngs.
11. Many conda environments are divergent across rules, especially for R plotting. However, it seems that some of them can be unified. Please unify them if possible.