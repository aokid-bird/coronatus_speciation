#==================#
#### Slice bams ####
#==================#
# Outgroup bam files are huge containing whole genome information, which is not required for
# the anaylsis of MIG-seq. THerefore, the bam files will be sliced based on the sites obtained from
# the MIG-seq analyses

#===========================#
#### Libraries & Sources ####
#===========================#
# libraries
library(tidyverse)
# sources
source("code/dataprep.R")

#======================#
#### 1) Directories ####
#======================#
#====================#
##### 1-1) Local #####
resdir_bwa.local <- "res/bwa"
resdir_angsd_intersect.local1 <- "res/angsd_intersect_mindepthind1"
resdir_angsd_intersect.local2 <- "res/angsd_intersect"
resdir_bwa_sra.local <- "res/bwa_sra"
#======================#
##### 1-2) Cluster #####
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")
resdir_angsd_intersect.clust1 <- str_interp("${dir_clust}/${resdir_angsd_intersect.local1}")
resdir_angsd_intersect.clust2 <- str_interp("${dir_clust}/${resdir_angsd_intersect.local2}")
resdir_bwa_sra.clust <- str_interp("${dir_clust}/res/bwa_sra")

#=========================#
#### 2) Obtained sites ####
#=========================#
# Because the intersecting sites may depend on the data subset, bam slicing will be done based on the
# largest set of sites obtained from the MIG-seq anlayses.
tmp.df <-
  data.df %>% 
  filter(!remove) %>% 
  mutate(bamlist.fullpath = str_c(resdir_bwa.clust, "/", sample, ".bam")) 
# population
pop <- levels(tmp.df$group_reg)

# setmindepthind = 1
geno1 <- 
  lapply(pop, 
         function(x){
           str_interp("zcat ${resdir_angsd_intersect.local1}/${x}/gl_reg.geno.gz") %>% 
             system(intern = TRUE) %>% 
             str_split("\t", simplify = TRUE) %>% 
             as_tibble()
         })
# setmindepthind = 2
geno2 <- 
  lapply(pop, 
         function(x){
           str_interp("zcat ${resdir_angsd_intersect.local2}/${x}/gl_reg.geno.gz") %>% 
             system(intern = TRUE) %>% 
             str_split("\t", simplify = TRUE) %>% 
             as_tibble()
         })
# largest set of the obtained sites
sites <-
  c(geno1, geno2) %>% 
  lapply(., select, V1, V2) %>% 
  lapply(., rename, scaf = V1, pos = V2) %>% 
  lapply(., 
         function(x){
           mutate(x, scaf_pos = str_c(scaf, ":", pos)) %>% 
             select(scaf_pos)
         }) %>% 
  bind_rows() %>% 
  pull(scaf_pos) %>% 
  unique() %>% 
  str_split(":", simplify = TRUE) %>% 
  as_tibble()
# save
write_delim(sites, 
            file = str_interp("${resdir_bwa.local}/slicing_sites.txt"),
            col_names = FALSE)

#===============================================================#
#### 4) Slice database bam files based on the obtained sites ####
#===============================================================#
#===================#
##### 4-1) Data #####
# paths for target bam files
ref.df <- 
  ref.df %>% 
  mutate(path.local = str_c(resdir_bwa_sra.local, "/", sraid, ".bam"),
         path.clust = str_c(resdir_bwa_sra.clust, "/", sraid, ".bam"),
  )

bamlist_ref <-
  ref.df %>% 
  filter(version %in% c("full", "combined")) %>% 
  pull(path.clust)

# obtained sites
sites <- 
  read_delim(str_interp("${resdir_bwa.local}/slicing_sites.txt"),col_names = FALSE) %>% 
  rename(scaf=X1, pos=X2)

#===========================#
##### 4-2) bam indexing #####
cddir <- str_interp("cd ${resdir_bwa_sra.clust}")
for(i in 1:length(bamlist_ref)){
  arg <- str_interp("samtools index ${bamlist_ref[i]}")
  
  cluster.scripter(filename = str_interp("sam_index_${bamlist_ref[i]}"),
                   local.directory = resdir_bwa_sra.local,
                   argument = list(cddir, arg),
                   name = str_interp("samindex${bamlist_ref[i]}"),
                   memsz = 250, cpunum = 72, envn = "bioinfo",
                   cluster.directory = resdir_bwa_sra.clust)
  
}

#=======================================#
##### 4-3) slice database bamfiles ######
# Database sequences were sliced based on the obtained sites
sites.group <- 
  sites %>% 
  group_split(scaf)

sites.clust <- 
  sites.group %>% 
  lapply(pull, pos) %>% 
  lapply(as.integer) %>% 
  lapply(group_consecutive)

scaf.group <-
  sites.group %>% 
  lapply(pull, scaf) %>% 
  lapply(unique)

position <- 
  map2(.x = scaf.group, 
       .y = sites.clust, 
       .f = function(x=.x,y=.y){str_c(x, ":", y)}) %>%
  lapply(str_c, collapse = " ") %>% 
  str_c(collapse = " ")
##### function end #####

#=========================#
##### 4-4) Call & Run #####
for(i in 1:length(bamlist_ref)){
  
  bam <- bamlist_ref[i]
  sample <- str_remove(bam, resdir_bwa_sra.clust) %>% str_remove_all("^/|.bam")
  
  bam.input.clust <- bam
  bam.output.clust <- bam %>% str_replace(".bam", str_interp("_slice.bam"))
  
  # samtools to view, sort and index bamfiles based on the intersecting sites
  arg.slice.all <- 
    str_interp("samtools view -o ${bam.output.clust} -O BAM -@ 71 ${bam.input.clust} ${position}")
  arg.sort.all <-
    str_interp("samtools sort -O BAM -o ${bam.output.clust} ${bam.output.clust}")
  arg.index.all <-
    str_interp("samtools index ${bam.output.clust}")
  
  cluster.scripter(filename = str_interp("bamslice_${sample}"),
                   local.directory = resdir_bwa_sra.local,
                   argument = list(arg.slice.all, arg.sort.all, arg.index.all),
                   name = str_interp("bsl${sample}"),
                   memsz = 250, cpunum = 72, envn = "bioinfo",
                   cluster.directory = resdir_bwa_sra.clust)
}