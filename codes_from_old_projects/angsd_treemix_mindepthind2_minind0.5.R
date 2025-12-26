#=================================================#
#### ANGSD for treemix analysis with glactools ####
#=================================================#
# This script intends to conduct treemix
# 1. We first re-computed GL for the intersecting sites.
# 2. Then, we created treemix input files using glactools.
# 3. Finally, conducting the treemix analyses on the allele frequency data

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
source("~/treemix-1.13/src/plotting_funcs.R")
library(OptM)
ource("code/dataprep.R")

rename <- dplyr::rename
mutate <- dplyr::mutate

data.df
ref.df

#======================#
#### 1) Directories ####
#======================#
#====================#
##### 1-1) Local #####
resdir_bwa.local <- "res/bwa"
resdir_angsd_intersect.local <- "res/angsd_intersect_mindepthind2"
resdir_angsd_treemix.local <- "res/angsd_treemix_mindepthind2"
resdir_bwa_sra.local <- "res/bwa_sra"

mkdir(resdir_angsd_treemix.local)

#======================#
##### 1-2) Cluster #####
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")
resdir_angsd_intersect.clust <- str_interp("${dir_clust}/${resdir_angsd_intersect.local}")
resdir_angsd_treemix.clust <- str_interp("${dir_clust}/${resdir_angsd_treemix.local}")
resdir_bwa_sra.clust <- str_interp("${dir_clust}/res/bwa_sra")

#================#
#### 2) ANGSD ####
#================#
#======================#
##### 2-1) Bamlist #####
file.copy(from = str_interp("${resdir_angsd_intersect.local}/bamlist_breed"),
          to = str_interp("${resdir_angsd_treemix.local}/bamlist"))
file.copy(from = str_interp("${resdir_angsd_intersect.local}/bamlist_breed_wout"),
          to = str_interp("${resdir_angsd_treemix.local}/bamlist_wout"))

#========================================#
##### 2-3) Remove Atlantisia rogersi #####
# remove Atlantisia
rmid.ref <- ref.df %>% filter(species == "Atlanticia rogersi") %>% pull(sraid)
read_lines(str_interp("${resdir_angsd_treemix.local}/bamlist_wout")) %>% str_subset(pattern = rmid.ref, negate = TRUE) %>% write_lines(str_interp("${resdir_angsd_treemix.local}/bamlist_wout"))

#===========================#
##### 2-3) linked sites #####
file.copy(from = str_interp("${resdir_angsd_intersect.local}/intersect_breed.txt"),
          to = str_interp("${resdir_angsd_treemix.local}/intersect_breed.txt"))
file.copy(from = str_interp("${resdir_angsd_intersect.local}/intersect_breed.chr"),
          to = str_interp("${resdir_angsd_treemix.local}/intersect_breed.chr"))

#=========================#
##### 2-4) Bamlist.df #####
bamlist.df <-
  bind_rows(
    data.df %>%
      filter(!remove) %>%
      filter(group_reg != "kanto_w") %>%
      mutate(ref = "data.df") %>%
      select(sample, group_reg, ref) %>% rename(group = group_reg),
    ref.df %>%
      filter(species == "Laterallus jamaicensis") %>%
      mutate(ref = "ref.df") %>%
      select(sraid, species, ref) %>% rename(sample = sraid, group = species)
  ) %>%
  mutate(bampath = case_when(
    ref == "data.df" ~ str_c(resdir_bwa.clust, "/", sample, ".bam"),
    ref == "ref.df"  ~ str_c(resdir_bwa_sra.clust, "/", sample, "_slice.bam")
  )) %>%
  mutate(group.label = str_remove_all(group, "_.*")) %>%
  rowwise() %>%
  mutate(group.label = case_when(
    str_detect(group.label, " ") ~ str_split(group.label, " ", simplify = TRUE) %>% str_sub(1,2) %>% str_c(collapse=""),
    TRUE ~ group.label
  )) %>%
  ungroup()

#====================================#
###### 2-4) parameter conditions #####
subdir <- 
  expand.grid(category = c("breed"),
              outgrp   = c("_wout")) %>% 
  mutate(subdir = str_c(category, outgrp)) %>% 
  pull(subdir)

#===============================#
##### 2-2) Angsd Parameters #####
###### 1:: do setting ######
gl             <- 1 # -GL
domaf          <- 1 # -doMaf
domajorminor   <- 1 # -domajorminor
doglf          <- 4 # -doGLF text loglikelihood
dogeno         <- 4 # -doGeno genotype txt format
docounts       <- 1 # -doCounts
dopost         <- 1 # -doPost
dobcf          <- 1 # -dovcf

###### 2:: filtering ######
minq <- 20
minmapq <- 30
# this time, both variant/invariant sites

uniqueonly <- 1
skiptriallelic <- 1
baq <- 2
###### 3:: misc ######
nthreads <- 10

#=========================#
##### 2-4) ANGSD CALL #####

for(i in 1:length(subdir)){
  
  tmp.subdir <- subdir[i]
  category <- "breed"
  sites <- str_interp("${resdir_angsd_treemix.clust}/intersect_breed.txt")
  chrs  <- str_interp("${resdir_angsd_treemix.clust}/intersect_breed.chr")
  outgrp <- ifelse(str_detect(tmp.subdir, "_wout"), "_wout", "")
  name.bamlist <-str_c("bamlist", outgrp)
  
  tmp.bamlist <- read_delim(str_interp("${resdir_angsd_treemix.local}/${name.bamlist}"), delim = " ", col_names = FALSE)
  
  dircreate <- str_interp("mkdir ${resdir_angsd_treemix.clust}/${tmp.subdir}")
  idx_rf <- str_interp("angsd sites index ${sites}")
  angsd <- str_c("angsd",
                 " -out ", str_interp("${resdir_angsd_treemix.clust}/${tmp.subdir}/gl_treemix"),
                 " -b ", str_interp("${resdir_angsd_treemix.clust}/${name.bamlist}"),
                 " -ref ", refpath_clust,
                 " -anc ", refpath_clust,
                 " -sites ", sites,
                 " -rf ", chrs,
                 " -GL ", gl,
                 #### DO parameters
                 " -doMaf ", domaf,
                 " -domajorminor ", domajorminor,
                 " -doGLF ", doglf,
                 " -doGeno ", dogeno,
                 " -doPost ", dopost,
                 " -doCounts ", docounts,
                 " -doBcf ", dobcf,
                 ### filtering parameters
                 " -uniqueOnly ", uniqueonly,
                 " -minMapQ ", minmapq,
                 " -minQ ", minq,
                 " -skipTriallelic ", skiptriallelic,
                 " -baq ", baq,
                 " -nthreads ", nthreads,
                 " 2> ", str_interp("${resdir_angsd_treemix.clust}/${tmp.subdir}/log_gl_treemix.oe")
  )
  
  cluster.scripter(filename = str_interp("angsd_treemix_${tmp.subdir}"),
                   local.directory = resdir_angsd_treemix.local,
                   argument = list(dircreate, idx_rf, angsd),
                   name = str_interp("agsd_trmx${tmp.subdir}"),
                   memsz = 50, cpunum = 10, envn = "angsd-0.94",
                   cluster.directory = resdir_angsd_treemix.clust)
}

#=========================#
####　3) Treemix input ####
#=========================#
subdir <- 
  expand.grid(category = c("breed"),
              outgrp   = c("_wout"),
              pop = c("pop", "indiv")) %>% 
  mutate(subdir = str_c(category, outgrp)) %>% 
  select(subdir, pop)

for(k in 1:length(unique(subdir$subdir))){
  
  tmp.subdir <- unique(subdir$subdir)[k]
  tmp.outg <- ifelse(str_detect(tmp.subdir, "wout"), "", "_wout")
  
  #=====================================#
  ##### 3-1) glactools to crate acf #####
  ###### 1:: Parameters ######
  fai <- str_replace(refpath_local, "fasta", "fasta.fai") %>% str_replace("..", "/mnt/y/Studies")
  header <- "#!/bin/bash"
  bcf <- list.files(str_interp("${resdir_angsd_treemix.local}/${tmp.subdir}"), pattern = ".bcf")
  vcf <- str_replace(bcf, "bcf", "vcf")
  vcf.simple <- str_replace(vcf, ".vcf", "_simple.vcf")
  acf <- str_replace(vcf, "vcf", "acf.gz")
  # 
  # ###### 2:: Call ######
  # 1 bcf2vcf via bcftools #
  call.bcf2vcf <- str_interp("bcftools convert -O v -o ${vcf} ${bcf}")
  # 2 remove PL/GP fields (only retainig GL for likelihoods) #
  call.rmpl <- str_interp("bcftools annotate -x FORMAT/GL ${vcf} > ${vcf.simple}")
  # 3 vcf2glf via glactools #
  call.vcf2acf <- str_interp("glactools vcfm2acf --fai ${fai} ${vcf.simple} > ${acf}")
  # 4 index acf
  call.index <- str_interp("glactools index ${acf}")
  
  ###### 3:: Write Bash #####
  list(header, call.bcf2vcf, call.rmpl, call.vcf2acf, call.index) %>%
    write_lines(str_interp("${resdir_angsd_treemix.local}/${tmp.subdir}/gen.acf.sh"))
  
  # do the following in the shell (in the ${resdir_angsd_treemix.local} directory)
  # bash gen.acf.sh
  
  for(t in 1:length(unique(subdir$pop))){
    tmp.pop <- as.character(unique(subdir$pop)[t])
    resdir_treemix <- str_interp("${resdir_angsd_treemix.local}/${tmp.subdir}/treemix_${tmp.pop}")
    mkdir(resdir_treemix)
    
    #========================#
    ##### 3-3) glactools #####
    pop <-
      str_interp("glactools view -p ${resdir_angsd_treemix.local}/${tmp.subdir}/${acf}") %>% system(intern = TRUE) %>%
      tibble(bampath = .)
    ###### 1:: rename outgroups ######
    rename.df <-
      bamlist.df %>%
      left_join(pop,., by="bampath") %>% 
      filter(!is.na(group.label))
    
    if(tmp.pop == "indiv"){
      rename.df %<>% 
        mutate(group.label = case_when(
          group == "Laterallus jamaicensis" ~ "Laja",
          TRUE ~ sample))
    }

    oldname <- str_c(rename.df$bampath, collapse = ",") %>% str_c("\"", ., "\"")
    newname <- str_c(rename.df$sample, collapse = ",") %>% str_c("\"", ., "\"")
    
    str_interp("glactools rename ${resdir_angsd_treemix.local}/${tmp.subdir}/${acf} ${oldname} ${newname} > ${resdir_treemix}/gl_treemix_tmp1.acf.gz") %>% system
    
    ##### 2:: merge populations ######
    merge.df <-
      rename.df %>% 
      group_by(group.label) %>% 
      group_split
    
    merge.call <-
      merge.df %>%
      lapply(., function(x){
        from <- pull(x, sample) %>%
          str_c(collapse = ",") %>%
          str_c("\"", ., "\"")
        to <- pull(x, group.label) %>%
          unique() %>%
          str_c("\"", ., "\"")
        str_c(from, " ", to)
      }) %>%
      str_c(collapse = " ")
    
    str_interp("glactools meld ${resdir_treemix}/gl_treemix_tmp1.acf.gz ${merge.call} > ${resdir_treemix}/gl_treemix_tmp2.acf.gz") %>% system
    
    if(str_detect(tmp.subdir, "wout")){
      ###### 3:: set root/anc populations ######
      str_interp("glactools usepopsrootanc ${resdir_treemix}/gl_treemix_tmp2.acf.gz Laja Laja > ${resdir_treemix}/gl_treemix_tmp3.acf.gz") %>%
        system
      str_interp("glactools replaceanc ${resdir_treemix}/gl_treemix_tmp2.acf.gz ${resdir_treemix}/gl_treemix_tmp3.acf.gz > ${resdir_treemix}/input.acf.gz") %>%
        system
    }else{
      system(str_interp("cp \"${resdir_treemix}/gl_treemix_tmp2.acf.gz\" \"${resdir_treemix}/input.acf.gz\""))
    }
    
    #=============================#
    ##### 3-4) Treemix input ######
    # treemix input with excluding population specific private alleles
    input.treemix <- str_interp("${resdir_treemix}/input.treemix")
    str_interp("glactools acf2treemix --noroot --noprivate ${resdir_treemix}/input.acf.gz > ${input.treemix} 2> ${resdir_treemix}/acf2treemix.log.oe") %>% system
    
    # version check
    str_interp("glactools view -P ${resdir_treemix}/input.acf.gz") %>% system
    # no version information
    
    #========================#
    ##### 3-5) Block Size ####
    #========================#
    resdir_lddecay <- str_interp("${resdir_treemix}/lddecay")
    mkdir(resdir_lddecay)
    blocksize.call <- str_interp("R --slave --vanilla --args ${resdir_angsd_global.local} ${tmp.subdir} ${resdir_treemix} < functions/lddecay_blocksizefinder.R 2> ${resdir_lddecay}/log.oe")
    system(blocksize.call)
    blocksize <- 
      read_lines(str_interp("${resdir_lddecay}/block_size.txt")) %>% 
      as.integer()
    
    #======================#
    #### 4) Treemix Run ####
    #======================#
    # The following block size parameter was defined in the R script "lddecay.R"
    # gunzip treemix input
    system(str_interp("gzip -c ${input.treemix} > ${input.treemix}.gz"))
    
    # loop treemix runs from migration edges = 0 to 6
    for(i in 1:6){
      for(o in 1:10){
        # number of edges to test
        edge = i-1
        # parameters
        numk <- edge # number of maximum migration edges
        blockk <- blocksize # block size
        if(str_detect(tmp.subdir, "wout")){
          outgroup <- "Laja"
          to_add_out <- str_interp(" -root ${outgroup}")
        }else{
          to_add_out <- ""
        }
        
        pathP <- "/home/aokid/miniconda3/envs/bioinfo/bin/consense"
        input <- str_interp("${resdir_treemix}/input.treemix.gz")
        outname <- str_interp("${resdir_treemix}/treemix_e${edge}_o${o}")
        logname <- str_interp("${resdir_treemix}/treemix_e${edge}_o${o}_log")
        # calls
        call.treemix <- str_interp("treemix -i ${input} -m ${numk} -o ${outname}${to_add_out} -bootstrap -k ${blockk} -se > ${logname}")
        
        system(call.treemix)
      }
    }
    
    #===========================#
    #### 5) Model evaluation ####
    #===========================#
    # evaluate the model with using delta m (Evanno method)
    optm <- OptM::optM(resdir_treemix, tsv = str_interp("${resdir_treemix}/optm.tsv"))
    OptM::plot_optM(optm, pdf = str_interp("${resdir_treemix}/optm.pdf"))
    
    #================================#
    ##### 5-1) Variance Explained ####
    # prepare lists that store data
    var_breed <- vector("list", 6)
    var_breed <- rep(list(var_breed),6)
    
    # calculate total variance explained
    for(i in 1:6){
      for(o in 1:10){
        edge <- i-1
        var_breed[[i]][[o]] <- c(unlist(calcVarExplain(str_interp("${resdir_treemix}/treemix_e${edge}_o${o}"))), edge = edge, o = o)
      }
    }
    
    # dataframe organization
    matvar <- var_breed %>%
      unlist(recursive = FALSE) %>%
      lapply(as.matrix) %>%
      lapply(t) %>%
      lapply(as_tibble) %>%
      bind_rows()
    
    matvar %>%
      mutate(edge = as.character(edge)) %>%
      group_by(edge) %>%
      summarise(mean = mean(VarExplain), sd = sd(VarExplain))
    
    
    #====================================================================#
    ##### 5-2) Log-likelihood and combine it with variance explained #####
    # likelihood & variance explained changes according to adding migration edges
    path <-
      list.files(resdir_treemix, pattern = "llik", full.names = TRUE, recursive = TRUE)
    
    df <-
      path %>%
      lapply(read_lines) %>%
      lapply(str_extract, "(?<=events: )[0-9.-]+(?= )") %>%
      unlist %>% matrix(ncol = 2, byrow=TRUE) %>%
      as_tibble
    
    df %<>% mutate(path = path) %>% mutate(edge = str_extract(path,"(?<=_e)[:digit:](?=_)") %>% as.numeric,
                                           rep  = str_extract(path, "[:digit:]+(?=.llik)")%>% as.numeric) %>%
      mutate(V1 = as.numeric(V1), V2 = as.numeric(V2))
    
    df %<>%
      left_join(matvar, by = c("edge" = "edge", "rep" = "o")) %>%
      mutate(diff = V2-V1)
    write_csv(df, file = str_interp("${resdir_treemix}/parameters.csv"))
    
    summary <-
      df %>%
      group_by(edge) %>%
      dplyr::summarise(mean = mean(V2, na.rm = TRUE), sd = sd(V2, na.rm = TRUE), se = sd/sqrt(length(.)),
                       mean.var = mean(VarExplain), sd.var = sd(VarExplain), se.var = sd.var/sqrt(length(.))) %>%
      ungroup %>%
      mutate(group = "rail")
    
    write_csv(summary, file = str_interp("${resdir_treemix}/result_summary.csv"))
    
    # pop information
    popordcol <-
      rename.df %>%
      mutate(group.label2 = str_remove_all(group, "_.*")) %>%
      mutate(group.label2 = factor(group.label2, levels = c("tomakomai", "kushiro", "aomori", "russia", "kanto", "Laja"))) %>% 
      arrange(group.label2) %>% 
      distinct(group.label)
    
    write_delim(popordcol, str_interp("${resdir_treemix}/popordcol"), delim = " ", col_names = FALSE)
  }
}
