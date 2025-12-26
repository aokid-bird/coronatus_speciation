#=======================#
#### QC and trimming ####
#=======================#
#
# This script conducts quality check (QC) of pre- and post- trimming
# Sequence trimming is based on the sequence quality and adapter sequences as well as PCR annealing sites

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
source("code/dataprep.R")

data.df

#==========================#
#### 1) Check raw fastq ####
#==========================#
#====================================#
##### 1-1) Directories & Sources #####
# primer/adapter reference file
ref.df <- read_csv("data/adapter_primer.csv")

#==================#
##### 1-2) Run #####
# use the adapter.check function in ngs_functions.R
fastq.check <-
  adapter.check(path = data.df$fullpath.R1[1],
                size = 1000,
                typelist = c("adapter", "ssr"), # these have to match with "type" variables in the ref
                ref = ref.df,
                option = "all") # a function is stored in "../functions/ngs_functions.R"

#======================#
##### 1-3) Results #####
# Visually check the results
id_adapter <- fastq.check %>% filter(type %in% c("pcrtail", "adapter")) %>% pull(id) %>% unique
id_ssr.rc <- fastq.check %>% filter(type == "ssr" & direction == "reverse complement") %>% pull(id) %>% unique

###### 1:: General cases with forward primers on the 5' side ######
fastq.check %>% filter(!((id %in% id_adapter) | (id %in% id_ssr.rc))) %>% nrow()
id_forward <- fastq.check %>% filter(!((id %in% id_adapter) | (id %in% id_ssr.rc))) %>% pull(id) %>% unique()

fastq.check %>% filter(id %in% id_forward)

# 996/1000 did not include full primer/adapter sequences

###### 2:: with adapters + ssr ######
fastq.check %>% filter(id %in% id_ssr.rc)

# 1/1000 included reverse complement adapter
# 3/1000 included reverse complement primers

#==========================#
#### 2) Pre-trimming QC ####
#==========================#
#====================================#
##### 2-1) Directories & Sources #####
resdir_fastqc.local <- "res/fastqc_pre"
resdir_fastqc.clust <- str_interp("${dir_clust}/${resdir_fastqc.local}")
if(!dir.exists(resdir_fastqc.local))dir.create(resdir_fastqc.local)

#=========================#
##### 2-2) Parameters #####
# contaminant files (MIG-seq PCR primers)
contaminants <- str_interp("data/primer.tsv")
# adapter files (ILLUMINA adapters)
adapters <- str_interp("data/adapter.tsv")

#=================#
##### 2-3) Run ####
# run fastqc
fastqc <- str_interp("fastqc -t 10 -o ./${resdir_fastqc.local} ")

for(i in 1:nrow(data.df)){
  df <- data.df[i,]
  sample <- df$sample
  sample.short <- str_extract(sample, ".*(?=_)")
  path <- df %>% select(fullpath.R1, fullpath.R2) %>% as.character

  for(k in 1:length(path)){

    fastqc.call <- str_interp("${fastqc}${path[k]} --contaminants ${contaminants} --adapters ${adapters}")

    system(fastqc.call)
  }

}

# run multiqc
multiqc <- str_interp("multiqc ${resdir_fastqc.local} --outdir ${resdir_fastqc.local}")
system(multiqc)

#======================#
#### 3) Trimmomatic ####
#======================#
#==========================#
##### 3-1) Directories #####
# local
resdir_trim.local <- "res/trimmomatic"

# cluster
resdir_trim.clust <- str_interp("${dir_clust}/${resdir_trim.local}")
datadir_adapter.clust <- str_interp("${dir_clust}/data")

# create a directory
if(!dir.exists(resdir_trim.local))dir.create(resdir_trim.local)

#=========================#
##### 3-2) Parameters #####
pe = "PE"
threads = 72
phred = 33
sw = "4:15"
# headcrop = 21
leading = 30
trailing = 30
minlen = 40
illuminaclip = "ILLUMINACLIP:$adapt/adapter_primer.fa:2:30:7"
#30 = palindromeClipThreshold, requiring 50 base matches between R1 and R2 to clip adapters
#7 = simpleClipThreshold, how accurate adapter sequence matches to the sequence. Low of the recommended (lower = more likely to match)

#==================#
##### 3-3) Run #####
for(i in 1:nrow(data.df)){
  x <- data.df[i,]

  # parameters
  prefix <- x$sample
  prefix.short <- str_extract(prefix, ".*(?=_)")
  path.R1 <- x$fullpath.R1_clust
  path.R2 <- x$fullpath.R2_clust

  # run call
  trimmomatic <- str_c("trimmomatic ", pe,
                       " -threads ", threads,
                       " -phred", phred,
                       " -trimlog ", resdir_trim.clust, "/", str_c(prefix, "_log.txt"),
                       " ", path.R1,
                       " ", path.R2,
                       " ", resdir_trim.clust, "/", prefix, "_pair_R1.fastq.gz",
                       " ", resdir_trim.clust, "/", prefix, "_unpair_R1.fastq.gz",
                       " ", resdir_trim.clust, "/", prefix, "_pair_R2.fastq.gz",
                       " ", resdir_trim.clust, "/", prefix, "_unpair_R2.fastq.gz",
                       # " HEADCROP:", headcrop, # try first with no headcrop and see how primers are removed by the adapter file
                       " ", illuminaclip,
                       " SLIDINGWINDOW:", sw,
                       " LEADING:", leading,
                       " TRAILING:", trailing,
                       " MINLEN:", minlen)
  # cat(trimmomatic)
  # system(trimmomatic)
  adapter <- str_interp("adapt=\"${datadir_adapter.clust}\"")

  # scripts for analysis on the cluster
  cluster.scripter(filename = str_interp("trim_${prefix}"),
                   local.directory = resdir_trim.local,
                   argument = list(adapter, trimmomatic),
                   name = str_interp("tr${prefix.short}"),
                   memsz = 250, cpunum = 72,
                   cluster.directory = resdir_trim.clust)
}

#======================#
##### 3-4) Cleanup #####
# create directories
# shell scripts
if(!dir.exists(str_interp("${resdir_trim.local}/shellscript")))dir.create(str_interp("${resdir_trim.local}/shellscript"))
# out/error files
if(!dir.exists(str_interp("${resdir_trim.local}/oe")))dir.create(str_interp("${resdir_trim.local}/oe"))
# move files - scripts
list.files(str_interp("${resdir_trim.local}"), pattern = "(que|run)_.*.sh", full.names = TRUE) %>%
  str_c(collapse = " ") %>%
  str_c("mv ", ., str_interp(" ${resdir_trim.local}/shellscript")) %>%
  system()

# move files - oe
list.files(str_interp("${resdir_trim.local}"), pattern = "\\.(o|e)[0-9]+$", full.names = TRUE) %>%
  str_c(collapse = " ") %>%
  str_c("mv ", ., str_interp(" ${resdir_trim.local}/oe")) %>%
  system()

#===========================#
#### 4) Post-trimming QC ####
#===========================#
#====================================#
##### 4-1) Directories & Sources #####
# local
resdir_fastqc.local <- "res/fastqc_post"
if(!dir.exists(resdir_fastqc.local))dir.create(resdir_fastqc.local)

#=========================#
##### 4-2) Parameters #####
# contaminant files (MIG-seq PCR primers)
contaminants <- str_interp("data/primer.tsv")
# adapter files (DNBSEQ adapters)
adapters <- str_interp("data/adapter.tsv")

#==================#
##### 4-3) Run #####
fastqc <- str_interp("fastqc -t 10 -o ./${resdir_fastqc.local} ")

for(i in 1:nrow(data.df)){
  df <- data.df[i,]
  sample <- df$sample
  sample.short <- str_extract(sample, ".*(?=_)")

  path <- df %>%
    mutate(R1_pair = str_c(resdir_trim.local, "/", sample, "_pair_R1.fastq.gz"),
           R2_pair = str_c(resdir_trim.local, "/", sample, "_pair_R2.fastq.gz")) %>%
    select(ends_with("pair")) %>%
    as.character

  for(k in 1:length(path)){

    fastqc.call <- str_interp("${fastqc}${path[k]} --contaminants ${contaminants} --adapters ${adapters}")

    system(fastqc.call)
  }

}

# run multiqc
multiqc <- str_interp("multiqc ${resdir_fastqc.local} --outdir ${resdir_fastqc.local}")
system(multiqc)
