#=======================#
#### ANGSD ABBA-BABA ####
#=======================#
# This script intends to calculate D-statistics using Abbababa2 program of ANGSD directly from GL.

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
library(BITEV2)
source("../functions/ngs_functions.R")
source("../functions/cluster_functions.R")
source("../functions/general_functions.R")
source("code/dataprep.R")

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
resdir_angsd_abbababa.local <- "res/angsd_abbababa_mindepthind2"
resdir_bwa_sra.local <- "res/bwa_sra"

mkdir(resdir_angsd_abbababa.local)

#======================#
##### 1-2) Cluster #####
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")
resdir_bwa_sra.clust <- str_interp("${dir_clust}/${resdir_bwa_sra.local}")
resdir_angsd_abbababa.clust <- str_interp("${dir_clust}/${resdir_angsd_abbababa.local}")

#================#
#### 2) ANGSD ####
#================#
#=========================#
##### 2-1) Parameters #####
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
  ungroup() %>% 
  mutate(group = factor(group, levels = c("tomakomai_b", "kushiro_b", "russia_b", "aomori_b", "Laterallus jamaicensis"))) %>% 
  arrange(group)

name.bamlist <- "bamlist"
bamlist <- 
  bamlist.df %>% 
  pull(bampath)

write_lines(bamlist, str_interp("${resdir_angsd_abbababa.local}/${name.bamlist}"))

#===========================#
##### 2-2) linked sites #####
# intersect sites from the intersecting run directory
file.copy(from = str_interp("${resdir_angsd_intersect.local}/intersect_breed.txt"),
          to = str_interp("${resdir_angsd_abbababa.local}/intersect_breed.txt"))
file.copy(from = str_interp("${resdir_angsd_intersect.local}/intersect_breed.txt.bin"),
          to = str_interp("${resdir_angsd_abbababa.local}/intersect_breed.txt.bin"))
file.copy(from = str_interp("${resdir_angsd_intersect.local}/intersect_breed.txt.idx"),
          to = str_interp("${resdir_angsd_abbababa.local}/intersect_breed.txt.idx"))
file.copy(from = str_interp("${resdir_angsd_intersect.local}/intersect_breed.chr"),
          to = str_interp("${resdir_angsd_abbababa.local}/intersect_breed.chr"))

#=======================#
##### 2-3) SizeFile #####
sizefile <- 
  bamlist.df %>% 
  group_by(group) %>% 
  summarize(n = n()) %>% 
  select(n)
write_delim(sizefile, str_interp("${resdir_angsd_abbababa.local}/sizeFile"), col_names = FALSE)

#=======================#
##### 2-4) popNames #####
sample.name <- 
  bamlist.df %>% 
  select(group.label) %>%
  mutate(group.label = case_when(
    group.label == "Laja" ~ "Laterallus_jamaicensis",
    TRUE ~ str_to_sentence(group.label)
  )) %>% 
  pull(group.label) %>% 
  unique

sample.name %>% write_lines(str_interp("${resdir_angsd_abbababa.local}/popNames"))

#==================================#
##### 2-5) Run angsd abbababa2 #####
###### 1:: Do parameters ######
docounts         <- 1 # -doCounts
doabbababa2          <- 1 # -doabbababa2

###### 2:: filtering ######
minq <- 20
minmapq <- 30
uniqueonly <- 1
skiptriallelic <- 1
baq <- 2

###### 3:: misc ######
uselast <- 1
nthreads <- 8
sites <- str_interp("${resdir_angsd_abbababa.clust}/intersect_breed.txt")
chrs  <- str_interp("${resdir_angsd_abbababa.clust}/intersect_breed.chr")

angsd <- str_c("angsd",
               " -doAbbababa2 ", doabbababa2,
               " -out ", str_interp("${resdir_angsd_abbababa.clust}/bam.Angsd"),
               " -b ", str_interp("${resdir_angsd_abbababa.clust}/${name.bamlist}"),
               " -sizeFile ", str_interp("${resdir_angsd_abbababa.clust}/sizeFile"),
               " -ref ", refpath_clust,
               " -anc ", refpath_clust,
               " -sites ", sites,
               " -rf ", chrs,
               " -useLast ", uselast,
               #### DO parameters
               " -doCounts ", docounts,
               ### filtering parameters
               " -uniqueOnly ", uniqueonly,
               " -minMapQ ", minmapq,
               " -minQ ", minq,
               # " -skipTriallelic ", skiptriallelic, # not accepted by abbababa
               " -baq ", baq,
               " -nthreads ", nthreads,
               " 2> ", str_interp("${resdir_angsd_abbababa.clust}/abbababa.oe")
)

cluster.scripter(filename = str_interp("angsd_abbababa"),
                 local.directory = resdir_angsd_abbababa.local,
                 argument = list(angsd),
                 name = str_interp("agsdabba"),
                 memsz = 50, cpunum = 10, envn = "angsd-0.94",
                 cluster.directory = resdir_angsd_abbababa.clust)

#=================#
#### 3) D-stat ####
#=================#
dstat <- str_interp("Rscript code/estAvgError.R angsdFile=\"${resdir_angsd_abbababa.local}/bam.Angsd\" out=\"${resdir_angsd_abbababa.local}/result\" sizeFile=\"${resdir_angsd_abbababa.local}/sizeFile\" nameFile=\"${resdir_angsd_abbababa.local}/popNames\"")

system(dstat)

result.Observed <- 
  read_table(str_interp("${resdir_angsd_abbababa.local}/result.Observed.txt")) %>% 
  mutate(file = "Uncorrected")
result.TransRem <- 
  read_table(str_interp("${resdir_angsd_abbababa.local}/result.TransRem.txt")) %>% 
  mutate(file = "Uncorrected Transition Removed")

npair <- nrow(result.Observed)

df.dstat <- 
  bind_rows(result.Observed, result.TransRem) %>% 
  mutate(sd = sqrt(`V(JK-D)`),
         adjsd = 3*sd,
         adjpvalue = pvalue*npair,
         label = case_when(
           adjpvalue < 0.05 & adjpvalue > 0.01 ~ "*",
           adjpvalue < 0.01 & adjpvalue > 0.001 ~ "**",
           adjpvalue < 0.001  ~ "***",
           TRUE ~ ""
         )) %>% 
  mutate(pair = str_c("D{", H1,",", H2, ":", H3, "}"))

write_csv(df.dstat, str_interp("${resdir_angsd_abbababa.local}/df.dstat.csv"))
