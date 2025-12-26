#=====================#
#### ANGSD for SFS ####
#=====================#
# This script intends to conduct SFS estimation using ANGSD for downstream demographic analyses
# Based on the PCA, ADMIXTURE, RAxML, NeighborNet, and SNAPP analyses, we grouped Aomori and Russia together as the "coontinental" group, and Tomakomai and Kushiro as the "island" group.
# 1. We first calculated genotype likelihood on the intersecting sites
# 2. Using realSFS of ANGSD,
#   Folded SFS
# were calculated for
#   linked all (both invariant/variant) sites

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
source("../functions/ngs_functions.R")
source("../functions/cluster_functions.R")
source("../functions/general_functions.R")
source("code/dataprep.R")

data.df
ref.df

#======================#
#### 1) Directories ####
#======================#

##### 1-1) Local #####
resdir_bwa.local <- "res/bwa"
resdir_angsd_intersect.local <- "res/angsd_intersect_mindepthind2"
resdir_angsd_sfs.local <- "res/angsd_sfs_mindepthind2"
resdir_bwa_sra.local <- "res/bwa_sra"

mkdir(resdir_angsd_sfs.local)

##### 1-2) Cluster #####
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")
resdir_angsd_intersect.clust <- str_interp("${dir_clust}/${resdir_angsd_intersect.local}")
resdir_angsd_sfs.clust <- str_interp("${dir_clust}/${resdir_angsd_sfs.local}")
resdir_bwa_sra.clust <- str_interp("${dir_clust}/res/bwa_sra")

#=======================================#
#### 2) ANGSD for intersecting sites ####
#=======================================#
##### 2-1) Parameters #####
###### 1:: group splitting ######
# bamlist.df <-
#   data.df %>%
#   filter(!remove) %>%
#   filter(group_reg != "kanto_w") %>%
#   mutate(ref = "data.df") %>%
#   select(sample, group_2, ref) %>%
#   rename(group = group_2) %>%
#   mutate(bampath = str_c(resdir_bwa.clust, "/", sample, ".bam"))

# ###### 2:: sub-directories #####
# subdir <-
#   bamlist.df %>%
#   pull(group) %>%
#   unique() %>%
#   as.character()

# for(i in 1:length(subdir)){
#   dir2mk <- str_interp("${resdir_angsd_sfs.local}/${subdir[i]}")
#   mkdir(dir2mk)
# }

# ###### 3:: bamlist #####
# bamlist <-
#   bamlist.df %>%
#   group_split(group) %>%
#   lapply(select, bampath, group)
# name.bamlist <- "bamlist"

# for(i in 1:length(bamlist)){
#   gr <- unique(bamlist[[i]]$group)
#   lst <- bamlist[[i]] %>% select(bampath)
#   write_delim(lst, str_interp("${resdir_angsd_sfs.local}/${gr}/${name.bamlist}"), delim = " ", col_names = FALSE)
# }

# ##### 2-2) Angsd Parameters #####
# ###### 1:: do setting ######
# gl           <- 1 # -GL
# domaf        <- 1 # -doMaf
# domajorminor <- 1 # -domajorminor
# doglf        <- 4 # -doGLF text loglikelihood
# dogeno       <- 4 # -doGeno genotype txt format
# docounts     <- 1 # -doCounts
# dopost       <- 1 # -doPost
# dobcf        <- 1 # -doBcf
# dosaf        <- 1 # -doSaf

# ###### 2:: filtering ######
# minq <- 20
# minmapq <- 30
# # no filter for SNP determination, snp_pval & minmaf & postcutoff

# uniqueonly <- 1
# skiptriallelic <- 1
# baq <- 2
# setmindepthind <- 2
# ###### 3:: misc ######
# nthreads <- 10

# # # region specific angsd
# for(i in 1:length(subdir)){
  
#   grp <- subdir[i]
#   tmp.bamlist <- bamlist.df %>% filter(group == grp)
#   minind <- floor(nrow(tmp.bamlist)*0.5) # at least half of the regional samples have the sites
  
#   angsd <- str_c("angsd",
#                  " -out ", str_interp("${resdir_angsd_sfs.clust}/${grp}/gl"),
#                  " -b ", str_interp("${resdir_angsd_sfs.clust}/${grp}/${name.bamlist}"),
#                  " -ref ", refpath_clust,
#                  " -anc ", refpath_clust,
#                  " -GL ", gl,
#                  " -doMaf ", domaf,
#                  " -domajorminor ", domajorminor,
#                  " -doGLF ", doglf,
#                  " -doGeno ", dogeno,
#                  " -doPost ", dopost,
#                  " -doCounts ", docounts,
#                  " -doBcf ", dobcf,
#                  " -dosaf ", dosaf,
#                  # no snp filtering
#                  " -uniqueOnly ", uniqueonly,
#                  " -setmindepthind ", setmindepthind,
#                  " -minMapQ ", minmapq,
#                  " -minQ ", minq,
#                  " -minInd ", minind,
#                  " -skipTriallelic ", skiptriallelic,
#                  " -baq ", baq,
#                  " -nthreads ", nthreads,
#                  " 2> ", str_interp("${resdir_angsd_sfs.clust}/${grp}/log_gl.oe")
#   )
  
#   cluster.scripter(filename = str_interp("angsd_sfs_${grp}"),
#                    local.directory = resdir_angsd_sfs.local,
#                    argument = list(angsd),
#                    name = str_interp("agsd_${grp}"),
#                    memsz = 50, cpunum = 10, envn = "angsd-0.94",
#                    cluster.directory = resdir_angsd_sfs.clust)
# }

# #=====================================#
# #### 3) Finding intersecting sites ####
# #=====================================#
# #===================#
# ##### 3-1) Data #####
# geno <-
#   lapply(subdir,
#          function(x){
#            str_interp("zcat ${resdir_angsd_sfs.local}/${x}/gl.geno.gz") %>%
#              system(intern = TRUE) %>%
#              str_split("\t", simplify = TRUE) %>%
#              as_tibble()
#          })

# # sites
# sites <-
#   geno %>%
#   lapply(., select, V1, V2) %>%
#   lapply(., rename, scaf = V1, pos = V2) %>%
#   lapply(., mutate, pos = as.integer(pos)) %>%
#   lapply(., arrange, scaf, pos)

# lapply(sites, nrow)
# #=========================================#
# ##### 3-2) Finding intersecting sites #####
# # intersecting sites among all the regions
# intersects <-
#   sites %>%
#   lapply(.,
#          function(x){
#            mutate(x, scaf_pos = str_c(scaf, ":", pos)) %>%
#              select(scaf_pos)
#          }) %>%
#   Reduce(intersect, .) %>%
#   separate(scaf_pos, c("scaf", "pos"), sep = ":") %>%
#   arrange(scaf, pos)
# nrow(intersects)
# # 130,368 including both variant/invariant sites
# # population-specific: insular 152,601, continental 168,143

# #=========================#
# ##### 3-3) Save files #####
# ###### 1:: sites #####
# write_delim(intersects,
#             file = str_interp("${resdir_angsd_sfs.local}/intersect.txt"),
#             col_names = FALSE)
# ###### 2:: scaf #####
# intersects %>%
#   select(scaf) %>%
#   distinct() %>%
#   write_delim(.,
#               file = str_interp("${resdir_angsd_sfs.local}/intersect.chr"),
#               col_names = FALSE)
# ###### 3:: index sites ######
# system(str_interp("angsd sites index ${resdir_angsd_sfs.local}/intersect.txt"))

# #=========================================#
# #### 4) Re-angsd on intersecting sites ####
# #=========================================#
# for(i in 1:length(subdir)){
  
#   grp <- subdir[i]
#   tmp.bamlist <- bamlist.df %>% filter(group == grp)
  
#   angsd <- str_c("angsd",
#                  " -out ", str_interp("${resdir_angsd_sfs.clust}/${grp}/gl_intersect"),
#                  " -b ", str_interp("${resdir_angsd_sfs.clust}/${grp}/${name.bamlist}"),
#                  " -ref ", refpath_clust,
#                  " -anc ", refpath_clust,
#                  " -sites ", str_interp("${resdir_angsd_sfs.clust}/intersect.txt"),
#                  " -rf ", str_interp("${resdir_angsd_sfs.clust}/intersect.chr"),
#                  " -GL ", gl,
#                  " -doMaf ", domaf,
#                  " -domajorminor ", domajorminor,
#                  " -doGLF ", doglf,
#                  " -doGeno ", dogeno,
#                  " -doPost ", dopost,
#                  " -doCounts ", docounts,
#                  " -doBcf ", dobcf,
#                  " -dosaf ", dosaf,
#                  # no snp filtering
#                  " -uniqueOnly ", uniqueonly,
#                  # no site filtering since already taken into account in the previous run
#                  " -minMapQ ", minmapq,
#                  " -minQ ", minq,
#                  " -baq ", baq,
#                  " -nthreads ", nthreads,
#                  " 2> ", str_interp("${resdir_angsd_sfs.clust}/${grp}/log_gl_intersect.oe")
#   )
  
#   cluster.scripter(filename = str_interp("angsd_sfs_${grp}_intersect"),
#                    local.directory = resdir_angsd_sfs.local,
#                    argument = list(angsd),
#                    name = str_interp("agsd${grp}int"),
#                    memsz = 50, cpunum = 10, envn = "angsd-0.94",
#                    cluster.directory = resdir_angsd_sfs.clust)
# }


#=========================#
#### 5) Estimating SFS ####
#=========================#
subdir.sfs <-
  expand.grid(pop = subdir,
              link = c("linked"),
              type = c("fold"))

for(i in 1:nrow(subdir.sfs)){
  grp <- subdir.sfs$pop[i]
  link <- ifelse(subdir.sfs$link == "linked", "", "_unlinked")[i]
  typ <- subdir.sfs$type[i]
  
  # parameters
  maxiter <- 50000
  tole <- "1e-6"
  fold <- ifelse(typ == "unfold", 0, 1) # 0=unfold, 1=fold
  thread <- 50
  optimsfs <- str_c("realSFS ",
                    str_interp("${resdir_angsd_sfs.clust}/${grp}/gl_intersect${link}.saf.idx"),
                    " -maxIter ", maxiter,
                    " -tole ", tole,
                    " -P ", thread,
                    " -fold ", fold,
                    " > ", str_interp("${resdir_angsd_sfs.clust}/${grp}/gl_intersect${link}_${typ}.sfs")
  )
  tajima <- str_c("realSFS saf2theta ",
                  str_interp("${resdir_angsd_sfs.clust}/${grp}/gl_intersect${link}.saf.idx"),
                  " -outname ",  str_interp("${resdir_angsd_sfs.clust}/${grp}/gl_intersect${link}_${typ}"),
                  " -sfs ",  str_interp("${resdir_angsd_sfs.clust}/${grp}/gl_intersect${link}_${typ}.sfs"),
                  ifelse(fold == 1, " -fold 1", "")
  )
  theta <- str_c("thetaStat do_stat ",
                 str_interp("${resdir_angsd_sfs.clust}/${grp}/gl_intersect${link}_${typ}.thetas.idx")
  )
  cluster.scripter(filename = str_interp("realsfs_${grp}${link}_${typ}"),
                   local.directory = resdir_angsd_sfs.local,
                   argument = list(optimsfs, tajima, theta),
                   name = str_interp("sfs${grp}${link}${typ}"),
                   memsz = 50, cpunum = 10, envn = "angsd-0.94",
                   cluster.directory = resdir_angsd_sfs.clust)
}

#====================================#
#### 7) 2D SFS estimation and FST ####
#====================================#
# rearrange the order for 2dSFS analysis (fastsimcoal2)
subdir <- factor(subdir, levels = c("continental", "insular"))
subdir <- subdir[order(subdir)] %>% as.character

subdir.sfs2d <-
  combn(subdir, 2) %>%
  t() %>%
  as_tibble() %>%
  rename(pair1 = V1, pair2 = V2) %>%
  expand_grid(link = c("linked"), type = c("fold"))

for(i in 1:nrow(subdir.sfs2d)){
  # parameters
  tmp.row <- subdir.sfs2d[i,]
  tmp.pair1 <- tmp.row$pair1
  tmp.pair2 <- tmp.row$pair2
  link <-ifelse(tmp.row$link == "linked", "", "_unlinked")
  typ <- tmp.row$type
  
  # run parameters
  maxiter <- 50000
  tole <- "1e-6"
  fold <- ifelse(typ == "unfold", 0, 1) # 0=unfold, 1=fold
  thread <- 10
  
  # files
  saf1     <- str_interp("${resdir_angsd_sfs.clust}/${tmp.pair1}/gl_intersect${link}.saf.idx")
  saf2     <- str_interp("${resdir_angsd_sfs.clust}/${tmp.pair2}/gl_intersect${link}.saf.idx")
  out.sfs  <- str_interp("${resdir_angsd_sfs.clust}/${tmp.pair1}_${tmp.pair2}_intersect${link}_${typ}.ml")
  out.bin  <- str_interp("${resdir_angsd_sfs.clust}/${tmp.pair1}_${tmp.pair2}_intersect${link}_${typ}.bin")
  out.fst  <- str_interp("${resdir_angsd_sfs.clust}/${tmp.pair1}.${tmp.pair2}_intersect${link}_${typ}_Fst")
  out.sw   <- str_interp("${resdir_angsd_sfs.clust}/${tmp.pair1}.${tmp.pair2}_intersect${link}_${typ}_sw")
  
  # call
  sfs2d    <- str_interp("realSFS ${saf1} ${saf2} -maxIter ${maxiter} -tole ${tole} -fold ${fold} -P ${thread} > ${out.sfs}")
  fst      <- str_c(str_interp("realSFS fst index ${saf1} ${saf2} -sfs ${out.sfs} -fstout ${out.fst}"),
                    ifelse(fold == 1, " -fold 1", ""))
  fst.glob <- str_interp("realSFS fst stats ${out.fst}.fst.idx > ${out.fst}.out")
  fst.window <- str_interp("realSFS fst stats2 ${out.fst}.fst.idx -win 50000 -step 10000 > ${out.sw]")
  sfsbin <- str_interp("realSFS bins ${saf1} ${saf2} -P 8 -fold ${fold} > ${out.bin}")
  
  cluster.scripter(filename = str_interp("realsfs2d_${tmp.pair1}_${tmp.pair2}${link}_${typ}"),
                   local.directory = resdir_angsd_sfs.local,
                   argument = list(sfs2d, fst, fst.glob,sfsbin),
                   name = str_interp("sfs2d${tmp.pair1}${tmp.pair2}${link}${typ}"),
                   memsz = 50, cpunum = 10, envn = "angsd-0.94",
                   cluster.directory = resdir_angsd_sfs.clust)
}

#===============================#
#### 8) Clean up directories ####
#===============================#
# list up files
# que/run files
qrf <- list.files(str_interp("${resdir_angsd_sfs.local}"), pattern = "(que|run)_*", full.names = TRUE)
# out/error files
oe  <- list.files(str_interp("${resdir_angsd_sfs.local}"), pattern = "(*\\.e[0-9]*|*\\.o[0-9]*)$*", full.names = TRUE)

# create directories
dir.qrf <- str_interp("${resdir_angsd_sfs.local}/que_run")
dir.oe  <- str_interp("${resdir_angsd_sfs.local}/oe")
mkdir(dir.qrf)
mkdir(dir.oe)

# copy and remove listed files
file.copy(from = qrf, to = dir.qrf, overwrite = TRUE)
file.remove(qrf)

file.copy(from = oe, to = dir.oe, overwrite = TRUE)
file.remove(oe)
