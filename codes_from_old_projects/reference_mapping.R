#==================================#
#### Reference mapping ILLUMINA ####
#==================================##
# This script maps trimmed/filtered sequences from ILLUMINA MIG-seq 
# on the reference genome of
# Laterallus jamaicensis by using the bwa software

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
source("code/dataprep.R")
data.df

#======================#
#### 1) Directories ####
#======================#
#====================#
##### 1-1) Local #####
resdir_trim.local <- "res/trimmomatic"
resdir_bwa.local <- "res/bwa"

mkdir(resdir_bwa.local)
#======================#
##### 1-2) Cluster #####
resdir_trim.clust <- str_interp("${dir_clust}/${resdir_trim.local}")
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")

#==============#
#### 2) Run ####
#==============#
# temporary file for paths
tmp <-
  data.df %>%
  mutate(path.pairR1.clust   = str_c(resdir_trim.clust, "/", sample, "_pair_R1.fastq.gz"),
         path.pairR2.clust   = str_c(resdir_trim.clust, "/", sample, "_pair_R2.fastq.gz"),
         path.unpairR1.clust = str_c(resdir_trim.clust, "/", sample, "_unpair_R1.fastq.gz"),
         path.unpairR2.clust = str_c(resdir_trim.clust, "/", sample, "_unpair_R2.fastq.gz")) %>%
  select(sample, starts_with("path."))

# bwa run and merge bam files
for(i in 1:nrow(data.df)){

  df <- tmp[i,]

  # parameters
  logname <- df$sample
  logname.short <- str_extract(logname, ".*(?=_)")
  algorithm = "mem"
  thread = 72

  # outfiles
  pbam <- str_interp("${resdir_bwa.clust}/${logname}_p.bam")
  unpbam1 <- str_interp("${resdir_bwa.clust}/${logname}_R1.bam")
  unpbam2 <- str_interp("${resdir_bwa.clust}/${logname}_R2.bam")
  out <- str_interp("${resdir_bwa.clust}/${logname}.bam")

  # bwas
  bwap <-
    str_c(
      str_interp("bwa ${algorithm} -t ${thread} ${refpath_clust}"),
      str_interp("${df$path.pairR1.clust} ${df$path.pairR2.clust} |"),
      str_interp("samtools sort -@ ${thread} -o ${pbam}"),
      sep = " ")

  bwa1 <-
    str_c(
      str_interp("bwa ${algorithm} -t ${thread} ${refpath_clust}"),
      str_interp("${df$path.unpairR1.clust} |"),
      str_interp("samtools sort -@ ${thread} -o ${unpbam1}"),
      sep = " ")

  bwa2 <-
    str_c(
      str_interp("bwa ${algorithm} -t ${thread} ${refpath_clust}"),
      str_interp("${df$path.unpairR2.clust} |"),
      str_interp("samtools sort -@ ${thread} -o ${unpbam2}"),
      sep = " ")

  # bam merge, sort and index
  merge <- str_interp("samtools merge ${out} ${pbam} ${unpbam1} ${unpbam2}")
  sort <- str_interp("samtools sort -@ ${thread} -o ${out} ${out}")
  index <- str_interp("samtools index ${out}")

  cluster.scripter(filename = str_interp("bwa_${logname}"),
                   local.directory = resdir_bwa.local,
                   argument = list(bwap, bwa1, bwa2, merge, sort, index),
                   name = str_interp("bw${logname.short}"),
                   memsz = 16, cpunum = 72,
                   cluster.directory = resdir_bwa.clust)
}

#==================#
#### 3) Cleanup ####
#==================#
# create directories
# shell scripts
if(!dir.exists(str_interp("${resdir_bwa.local}/shellscript")))dir.create(str_interp("${resdir_bwa.local}/shellscript"))
# out/error files
if(!dir.exists(str_interp("${resdir_bwa.local}/oe")))dir.create(str_interp("${resdir_bwa.local}/oe"))
# move files - scripts
list.files(str_interp("${resdir_bwa.local}"), pattern = "(que|run)_.*.sh", full.names = TRUE) %>%
  str_c(collapse = " ") %>%
  str_c("mv ", ., str_interp(" ${resdir_bwa.local}/shellscript")) %>%
  system()

# move files - oe
list.files(str_interp("${resdir_bwa.local}"), pattern = "\\.(o|e)[0-9]+$", full.names = TRUE) %>%
  str_c(collapse = " ") %>%
  str_c("mv ", ., str_interp(" ${resdir_trim.local}/oe")) %>%
  system()
