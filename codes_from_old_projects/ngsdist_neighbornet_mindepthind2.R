#=============================#
#### ngsDist & NeighborNet ####
#=============================#
# This script intends to conduct distAngsd, whose output is used in NeighborNet analysis

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
library(tanggle)
library(ggtree)
source("code/dataprep.R")
source("../functions/ngs_functions.R")
source("../functions/cluster_functions.R")
source("../functions/general_functions.R")

data.df
ref.df

#======================#
#### 1) Directories ####
#======================#
##### 1-1) Local #####
resdir_dist.local <- "res/ngsdist_mindepthind2"
resdir_angsd.local <- "res/angsd_global_mindepthind2"
resdir_angsd_run.local <- str_interp("${resdir_angsd.local}/breed_wout")
resdir_bwa.local <- "res/bwa"
resdir_bwa_sra.local <- "res/bwa_sra"
mkdir(resdir_dist.local)
#####1-2) Cluster #####
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")
resdir_bwa_sra.clust <- str_interp("${dir_clust}/${resdir_bwa_sra.local}")
#==================#
#### 2) ngsDist ####
#==================#
##### 2-1) names #####
unlinkedbeagle <- str_interp("${resdir_angsd_run.local}/gl_global_unlinked.beagle.gz")
distout <- str_interp("${resdir_dist.local}/ngsdist")
posinfo <- str_interp("${resdir_dist.local}/posinfo")
label <- str_interp("${resdir_dist.local}/label")
# out for p-distance
distoutp <- str_interp("${distout}_p")
distlogoutp <- str_interp("${distout}_p.oe")
# out for JC69
distoutjc <- str_interp("${distout}_jc69")
distlogoutjc <- str_interp("${distout}_jc69.oe")

##### 2-2) Parameters #####
# bamlist
bamlist <- 
  read_lines(str_interp("${resdir_angsd.local}/bamlist_breed_wout")) %>% 
  tibble %>% 
  set_colnames("bamlist.fullpath") %>% 
  mutate(sample = str_remove_all(bamlist.fullpath, str_interp("${resdir_bwa.clust}/|.bam|_slice|${resdir_bwa_sra.clust}/"))) %>% 
  left_join(., data.df, by = "sample")
# label
label.line <-
  bamlist %>% 
  pull(sample) 
# beagle & sites
beagle <- 
  system(str_interp("zcat ${unlinkedbeagle}"), intern=TRUE) %>% 
  str_split("\t", simplify=TRUE) %>% as_tibble %>% slice(-1)
nsites <- nrow(beagle)
posinfo.df <- 
  beagle %>% 
  select(V1, V2, V3) %>% 
  separate(V1, c("chr", "site"), sep = "_")

write_lines(label.line, label)
write_tsv(posinfo.df, posinfo, col_names = FALSE)

##### 2-3) Runs #####
# scripts - p-distance
call.ngsdist <- str_interp("~/ngsTools/ngsDist/ngsDist --probs true --avg_nuc_dist --evol_model 0 --labels ${label} --geno ${unlinkedbeagle} --n_ind ${nrow(bamlist)} --n_sites ${nsites} --out ${distoutp} 2> ${distlogoutp}")
system(call.ngsdist)

# scripts - JC69
call.ngsdist <- str_interp("~/ngsTools/ngsDist/ngsDist --probs true --avg_nuc_dist --evol_model 1 --labels ${label} --pos ${posinfo} --geno ${unlinkedbeagle} --n_ind ${nrow(bamlist)} --n_sites ${nsites} --out ${distoutjc} 2> ${distlogoutjc}")
system(call.ngsdist)

#======================#
#### 3) Nexus input ####
#======================#
# Create a nexus input file from ngsdist result
# p-distance
mat <- 
  read_tsv(distoutp, skip = 2, col_names = FALSE) %>% 
  rename(ind = X1)
treeinput_p <- str_interp("${resdir_dist.local}/ngsdist_input_p.nexus")
dist2nexus(df=mat, dryrun = FALSE, dir = treeinput_p, labels = "yes", triangle = "both")

# JC69
mat <- 
  read_tsv(distoutjc, skip = 2, col_names = FALSE) %>% 
  rename(ind = X1) 
treeinput_jc <- str_interp("${resdir_dist.local}/ngsdist_input_jc.nexus")
dist2nexus(df=mat, dryrun = FALSE, dir = treeinput_jc, labels = "yes", triangle = "both")

#======================#
#### 4) NeighborNet ####
#======================#
# reconstruct neighbornet from the nexus input
# p-distance
dm_p <- phangorn::read.nexus.dist(treeinput_p)
# Nnet_p <- phangorn::neighborNet(dm_p)
# read the output of SplitsTree4
Nnet_p <- phangorn::read.nexus.networx(str_interp("${resdir_dist.local}/nnet_p.nexus"))
plot(Nnet_p)
# jc69
dm_jc <- phangorn::read.nexus.dist(treeinput_jc)
# Nnet_jc <- phangorn::neighborNet(dm_jc)
Nnet_jc <- phangorn::read.nexus.networx(str_interp("${resdir_dist.local}/nnet_jc.nexus"))
