#================================================================#
#### ANGSD for SNAPP analysis and SNAPP species tree estimate ####
#================================================================#
# This script intends to conduct analyses on SNAPP
# using the breeding dataset. ANGSD is conducted specifically on only samples that scored lowest missingness
# which is approx four samples for each population.
# 1. We first determined the samples to be used for each population.
# 2. Then, we conducted ANGSD on these samples with minInd = max(sample size)
# 3. Slice the bamfile for outgroups
# 4. And finally call genotypes by using ANGSD and create a nexus file
# 5. and conduct SNAPP using the dataset

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
library(phangorn)
library(ggrepel)
library(treeio)
library(tidytree)
library(ggtree)
library(khroma)
source("../functions/ngs_functions.R")
source("../functions/cluster_functions.R")
source("../functions/general_functions.R")
source("code/dataprep.R")

data.df
ref.df

#======================#
#### 1) Directories ####
#======================#
#====================#
##### 1-1) Local #####
resdir_bwa.local <- "res/bwa"
resdir_angsd_intersect.local <- "res/angsd_intersect_mindepthind2"
resdir_angsd_global.local <- "res/angsd_global_rlrmv_mindepthind2"
resdir_angsd_snapp.local <- "res/angsd_snapp_mindepthind2"
resdir_bwa_sra.local <- "res/bwa_sra"

mkdir(resdir_angsd_snapp.local)

#======================#
##### 1-2) Cluster #####
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")
resdir_angsd_intersect.clust <- str_interp("${dir_clust}/${resdir_angsd_intersect.local}")
resdir_angsd_global.clust <- str_interp("${dir_clust}/${resdir_angsd_intersect.local}")
resdir_angsd_snapp.clust <- str_interp("${dir_clust}/${resdir_angsd_snapp.local}")
resdir_bwa_sra.clust <- str_interp("${dir_clust}/res/bwa_sra")

#===========================#
#### 2) Sample selection ####
#===========================#
gldir <- str_interp("${resdir_angsd_global.local}/breed")
bamlist.path <- str_interp("${resdir_angsd_global.local}/bamlist_breed")
geno.path <- str_interp("${gldir}/gl_global.geno.gz")
name <- read_lines(bamlist.path) %>% str_remove_all(str_interp("${resdir_bwa.clust}/|${resdir_bwa_sra.clust}/|.bam"))
geno.df <- 
  system(str_interp("zcat ${geno.path}"), intern = TRUE) %>% 
  str_split("\t", simplify = TRUE) %>% 
  .[,-ncol(.)] %>% 
  as_tibble() %>% 
  set_colnames(c("scaf", "pos", name)) %>% 
  mutate(pos = as.integer(pos)) %>% 
  dplyr::arrange(scaf, pos) %>% 
  unite(col = "scaf_pos", scaf, pos)
missing.df <- 
  geno.df %>% 
  select(!scaf_pos) %>% 
  mutate_all(str_count, "NN") %>% 
  apply(2, sum) %>% 
  tibble(sample = names(.), n.miss = .) %>% 
  mutate(prop.miss = n.miss/nrow(geno.df))
summary_df <- 
  left_join(missing.df, data.df, by = "sample") %>% 
  select(sample, group_reg, n.miss, prop.miss) 
# check samples
ggplot(data = summary_df, aes(x = group_reg, y = prop.miss))+
  geom_point(aes(color = log(prop.miss)))+
  geom_text_repel(aes(label = sample))+
  theme(legend.position = "none")+
  scale_color_viridis_c()
ggsave(str_interp("${resdir_angsd_snapp.local}/sample_prop.miss.png"))
write_csv(summary_df, str_interp("${resdir_angsd_snapp.local}/summary_df.csv"))

# select best 3-4 samples
samples_df <- 
  summary_df %>%
  group_by(group_reg) %>% 
  arrange(prop.miss) %>% 
  slice_head(n = 4) %>% 
  ungroup %>% 
  mutate(REMOVE = (prop.miss > 0.1))
# save
samples_df %>% 
  filter(!REMOVE) %>% 
  pull(sample) %>% 
  write_lines(str_interp("${resdir_angsd_snapp.local}/sample.txt"))

#================#
#### 3) ANGSD ####
#================#
#======================#
##### 3-1) Bamlist #####
file.copy(from = str_interp("${resdir_angsd_global.local}/bamlist_breed_wout"),
          to = str_interp("${resdir_angsd_snapp.local}/bamlist"),
          overwrite = TRUE)
# Retain the best samples 
retainedid <- read_lines(str_interp("${resdir_angsd_snapp.local}/sample.txt"))

# Laterallus as an outgroup
tmp <- 
  ref.df %>%
  filter(species == "Laterallus jamaicensis") %>%
  select(sraid) %>%
  mutate(bamlist.fullpath = str_c(dir_clust, "/", resdir_bwa_sra.local, "/", sraid, "_slice.bam")) %>%
  select(bamlist.fullpath) %>% 
  as.character()

# remove ones of pairs of related individuals 
read_lines(str_interp("${resdir_angsd_snapp.local}/bamlist")) %>% 
  str_subset(pattern = str_c(retainedid, collapse = "|"), negate = FALSE) %>% 
  c(., tmp) %>% 
  write_lines(str_interp("${resdir_angsd_snapp.local}/bamlist"))
bamlist <- 
  read_lines(str_interp("${resdir_angsd_snapp.local}/bamlist"))
#=========================#
##### 3-2) Sites/Scaf #####
# copy from the angsd_intersect
sites.original <- str_interp("#{resdir_angsd_snapp.local}/intersect_breed.txt")
scaf.original <- str_interp("#{resdir_angsd_snapp.local}/intersect_breed.chr")
sites.local <- str_interp("${resdir_angsd_snapp.local}/intersect_breed.txt")
scaf.local <- str_interp("${resdir_angsd_snapp.local}/intersect_breed.chr")
file.copy(sites.original, sites.local)
file.copy(sacf.original, scaf.local)
sites.clust  <- str_interp("${resdir_angsd_snapp.clust}/intersect_breed.txt")
scaf.clust   <- str_interp("${resdir_angsd_snapp.clust}/intersect_breed.chr")

#===============================#
##### 3-4) Angsd Parameters #####
###### 1:: do setting ######
gl           <- 1 # -GL samtools
domaf        <- 1 # -doMaf
domajorminor <- 1 # -domajorminor
doglf        <- 2 # -doGLF text beagle
dogeno       <- 4 # -doGeno genotype txt format
docounts     <- 1 # -doCounts
dopost       <- 1 # -doPost
dobcf        <- 1 # bcf

###### 2:: filtering ######
minq <- 20
minmapq <- 30
minind <- length(bamlist)
minmaf <- 0.1
snp_pval <- 1e-06
skiptriallelic <- 1

uniqueonly <- 1
baq <- 2
postcutoff <- 0.9 
# this time more confidence for snp calling

# mindepthind
# was not included here since the LD unlinked sites were as results of these 
# filters in the "global" analysis

###### 3:: misc ######
nthreads <- 10
##### 4:: angsd call #####
idx.sites <- str_interp("angsd sites index ${sites.clust}")
angsd <- str_c("angsd",
               " -out ", str_interp("${resdir_angsd_snapp.clust}/gl"),
               " -b ",   str_interp("${resdir_angsd_snapp.clust}/bamlist"),
               " -ref ", refpath_clust,
               " -anc ", refpath_clust,
               " -GL ", gl,
               " -sites ", sites.clust,
               " -rf ", scaf.clust,
               ### do option ###
               " -doMaf ", domaf,
               " -domajorminor ", domajorminor,
               " -doGLF ", doglf,
               " -doGeno ", dogeno,
               " -doPost ", dopost,
               " -doCounts ", docounts,
               " -doBcf ", dobcf,
               ### read/mapping quality filter option ###
               " -uniqueOnly ", uniqueonly,
               " -minMapQ ", minmapq,
               " -minQ ", minq,
               " -skipTriallelic ", skiptriallelic,
               " -baq ", baq,
               ### site filter option ###
               " -SNP_pval ", snp_pval,
               " -minInd ", minind,
               " -minMaf ", minmaf,
               ### snp calling option ###
               " -postCutoff ", postcutoff,
               # postcutoff included for calling snpgs
               " -nthreads ", nthreads,
               " 2> ", str_interp("${resdir_angsd_snapp.clust}/log_gl.oe")
)

cluster.scripter(filename = str_interp("angsd_snapp"),
                 local.directory = resdir_angsd_snapp.local,
                 argument = list(idx.sites, angsd),
                 name = str_interp("agsdsnapp"),
                 memsz = 50, cpunum = 12, envn = "angsd-0.94",
                 cluster.directory = resdir_angsd_snapp.clust)

#=====================#
#### 4) LD pruning ####
#=====================#
# Due to the reduction of the sample size, LD could be overestimated.
# Therefore, ngsLD results in angsd_global will be used to prune the sites
# obtained here

file.copy(str_interp("${gldir}/LD_unlinked.id"),
          str_interp("${resdir_angsd_snapp.local}/unlinkedsnp_pre.txt"))

unlinkedsites <- read_lines(str_interp("${resdir_angsd_snapp.local}/unlinkedsnp_pre.txt"))

# check how many sites to be retained
tmp <- 
  system(str_interp("zcat ${resdir_angsd_snapp.local}/gl.beagle.gz"), intern = TRUE) %>%
  str_split("\t", simplify = TRUE) %>% 
  as_tibble()

colnames(tmp) <- tmp[1,]
tmp <- tmp[-1,]

length(which(tmp$marker %in% str_replace(unlinkedsites, ":", "_")));nrow(tmp)
# 246/307 sites

unlinked <- 
  tmp$marker[tmp$marker %in% str_replace(unlinkedsites, ":", "_")] %>% 
  str_split("_", simplify = TRUE) %>% 
  as_tibble() %>% 
  magrittr::set_colnames(c("scaf", "pos")) %>% 
  arrange(scaf, pos)
sites <- str_interp("${resdir_angsd_snapp.local}/unlinkedsnp.txt")
write_tsv(unlinked, sites, col_names = FALSE)

#=====================#
#### 5) Nexus prep ####
#=====================#
#===============================================#
##### 5-1) Filter out linked sites from vcf #####
# file stem
gz <- str_interp("${resdir_angsd_snapp.local}/gl")

# population file
popref <-
  bind_rows(data.df %>% select(sample, group_reg) %>% mutate(group_reg = str_remove(group_reg, "_b|_w")) %>% set_colnames(c("sample", "pop")),
            ref.df %>% select(sraid, species) %>% mutate(species = str_replace(species, " ", "")) %>% set_colnames(c("sample", "pop")))

rename.df <-
  tibble(bamlist.fullpath = bamlist) %>%
  mutate(sample = str_remove(bamlist.fullpath, str_c(resdir_bwa.clust, "/|", resdir_bwa_sra.clust,"/")) %>% str_remove(".bam|_slice.*.bam")) %>%
  left_join(.,popref, by="sample") %>%
  group_by(pop) %>%
  mutate(n = row_number()) %>%
  ungroup %>%
  mutate(pop_n = str_c(pop, "_", n)) %>%
  select(pop_n)

sample <- str_interp("${resdir_angsd_snapp.local}/sample.txt")
write_delim(rename.df, sample,col_names = FALSE)

# create excess-heterozygous/linked-sites-filtered vcf file
bcf.hetfilt <- str_interp("bcftools +fill-tags ${gz}.bcf -Ou -- -t HWE | bcftools view -Oz -e'HWE<=0.05' -T ${sites} | bcftools reheader -s ${sample} > ${gz}_hetfilt.vcf.gz")
system(bcf.hetfilt)

# check if subsetting succeeded
n <- system(str_interp("bcftools query -f 'x' ${gz}_hetfilt.vcf.gz | wc -c"), intern=TRUE)
n == nrow(unlinked) # FALSE becayse hetfilt was done priori
n.hetfilt <- system(str_interp("bcftools query -f 'x' ${gz}_hetfilt.vcf.gz | wc -c"), intern=TRUE)
cat(str_interp("# unlinked sites with no heterozygosity excess = ${n.hetfilt}. # of all unlinked sites = ${nrow(unlinked)}"))

# convert vcf to phylip
system(str_interp("~/vcf2phylip/vcf2phylip.py -i ${gz}_hetfilt.vcf.gz --output-folder ${resdir_angsd_snapp.local} -o Laterallusjamaicensis_1 --nexus-binary"))

#=============================================#
#### 6) SNAPP analysis using snapp_prep.rb ####
#=============================================#
# snapp_prep.rb implements divergence time estimate by calibration
# Different prior for calibration time may be needed
# This script requires population/species definition file, constraint file (time calibration)

##### 6-1) File prep - pop file #####
pop.df <-
  rename.df %>%
  mutate(species = str_extract(pop_n, ".*(?=_)")) %>%
  rename(individual = pop_n) %>%
  select(species, individual)
write_tsv(pop.df, str_interp("${resdir_angsd_snapp.local}/pop.txt"))

##### 6-2) File prep - constraint file #####
# constraint[1,]: determined by the mitochondrial tree reconstructed in this study using molecular clock
# constraint[2,]: determined by the previous fossil calibrated trees (Garcia-R 2020 Fig2) to include the range of both Fig2A and Fig2B
constraint <-
  tibble(dist = c("lognormal(0,5.11,0.3)", "lognormal(0,8.95,0.15)"),
         type  = c("crown", "crown"),
         taxon = c("Laterallusjamaicensis,tomakomai,kushiro,russia,aomori", "Laterallusjamaicensis,tomakomai,kushiro,russia,aomori"))
write_tsv(constraint[1,], str_interp("${resdir_angsd_snapp.local}/constraint.txt"), col_names = FALSE)
write_tsv(constraint[2,], str_interp("${resdir_angsd_snapp.local}/constraint2.txt"), col_names = FALSE)

#======================================#
##### 6-3) Call for xml generation #####
call <- str_interp("ruby ~/snapp_prep.rb -p ${resdir_angsd_snapp.local}/gl_hetfilt.min4.phy -t ${resdir_angsd_snapp.local}/pop.txt -c ${resdir_angsd_snapp.local}/constraint.txt -l 500000 -o snapp_heightfix -x ${resdir_angsd_snapp.local}/snapp_heightfix.xml")
call2 <- str_interp("ruby ~/snapp_prep.rb -p ${resdir_angsd_snapp.local}/gl_hetfilt.min4.phy -t ${resdir_angsd_snapp.local}/pop.txt -c ${resdir_angsd_snapp.local}/constraint2.txt -l 500000 -o snapp_heightfix2 -x ${resdir_angsd_snapp.local}/snapp_heightfix2.xml")
system(call)
system(call2)

##### 6-4) SNAPP analysis on BEAST #####
# Used GUI BEAST2 v.2.7.3 for run

#==============================#
#### 7) SNAPP tree Analysis ####
#==============================#
# Analyze the tree output by snapp tree set analyzer app of BEAST and 
# manually save the output as a text file in the directory
##### out1 #####
out1 <- read_lines(str_interp("${resdir_angsd_snapp.local}/snapp_tree_set_analyzer_heightfix.txt"), skip = 3)
out1 <- out1[-(which(out1 == ""):length(out1))]
# percentage data
summary <- 
  str_extract(out1, ".*%") %>% 
  str_split(": ", simplify = TRUE) %>% 
  as_tibble %>% 
  magrittr::set_colnames(c("treeid", "percentage")) %>% 
  mutate(treeid = str_replace(treeid, " ", "") %>% str_to_lower(),
         percentage = str_remove(percentage, "%") %>% as.numeric)
out1 <- str_remove(out1, ".*% ")

# write out each as a single newick
for(i in 1:length(out1)){
  write_lines(out1[i], str_interp("${resdir_angsd_snapp.local}/hpd_heightfix_tree${i}.nwk"))
}
# save summary df
write_csv(summary, str_interp("${resdir_angsd_snapp.local}/hpd_heightfix_summary.csv"))

# treeanotation for 95%HPD trees
for(i in 1:length(out1)){
  system(str_interp("/home/aokidaisuke/miniconda3/envs/bioinfo/bin/treeannotator -heights median -burnin 50000 -target ${resdir_angsd_snapp.local}/hpd_heightfix_tree${i}.nwk ${resdir_angsd_snapp.local}/snapp_heightfix.trees ${resdir_angsd_snapp.local}/snapp_heightfix_annotated${i}.trees"))
}

# treeannotation for maximum clade credibility tree
system(str_interp("/home/aokidaisuke/miniconda3/envs/bioinfo/bin/treeannotator -heights median -burnin 50000 ${resdir_angsd_snapp.local}/snapp_heightfix.trees ${resdir_angsd_snapp.local}/snapp_heightfix_annotatedMCC.trees"))

levels <- c("tomakomai", "kushiro", "russia", "aomori", "outgroup")

tree.df <- vector("list", length(out1))
for(i in 1:length(out1)){
  tree.df[[i]] <- read.beast(str_interp("${resdir_angsd_snapp.local}/snapp_heightfix_annotated${i}.trees")) %>% as_tibble() 
  tree.df[[i]] <- 
    tree.df[[i]] %>% 
    mutate(
      percentage = summary$percentage[i],
      HPD_1 = map(height_0.95_HPD, .f = function(x)x[1]) %>% unlist %>% round(., 2),
      HPD_2 = map(height_0.95_HPD, .f = function(x)x[2]) %>% unlist %>% round(., 2),
      posterior.round = round(posterior, 2),
      height_median.round = signif(height_median, 3),
      posterior_95HPD = str_c(posterior.round, "\n", height_median.round, " [", HPD_1, ", ", HPD_2, "]")) %>% 
    mutate(label = if_else(str_detect(label, "Laterallus"), "outgroup", label) %>% factor(., levels = levels))
  tree.df[[i]] <- as.treedata(tree.df[[i]])
}

tree.mcc.df <- read.beast(str_interp("${resdir_angsd_snapp.local}/snapp_heightfix_annotatedMCC.trees")) 

colv <- c(as.character(color("bright")(6))[c(2,4,1,3)], "#4D4D4D")
names(colv) <- levels

p.tree <- 
  ggdensitree(tree.df,aes(alpha = percentage*.01, color = percentage), size = 2) +
  geom_tippoint(aes(fill = label), pch = 21, size = 5)+
  geom_tiplab(aes(label=str_to_title(label)), size = 5, nudge_x = .3)+
  geom_nodelab(aes(x=-height_median, label=posterior_95HPD),hjust = 1, vjust = -.5, size=5) +
  geom_range(range = 'height_0.95_HPD', color = "grey66", alpha = .6, size = 2)+
  coord_cartesian(clip="off")+
  scale_color_viridis_c("Percentage", option = "B", direction = -1)+
  scale_fill_manual(values = colv)+
  guides(fill = "none", alpha = "none")+
  theme(
    plot.margin = unit(c(1, 3, 1, 1), "cm") ,
    legend.position = "bottom"
  )

p.tree
ggsave(str_interp("figures/revise/Fig2/snapp_hightfix1.png"))
ggsave(str_interp("figures/revise/Fig2/snapp_hightfix1.eps"))

##### out2 #####
out2 <- read_lines(str_interp("${resdir_angsd_snapp.local}/snapp_tree_set_analyzer_heightfix2.txt"), skip = 3)
out2 <- out2[-(which(out2 == ""):length(out2))]
# percentage data
summary <- 
  str_extract(out2, ".*%") %>% 
  str_split(": ", simplify = TRUE) %>% 
  as_tibble %>% 
  magrittr::set_colnames(c("treeid", "percentage")) %>% 
  mutate(treeid = str_replace(treeid, " ", "") %>% str_to_lower(),
         percentage = str_remove(percentage, "%") %>% as.numeric)
out2 <- str_remove(out2, ".*% ")

# write out each as a single newick
for(i in 1:length(out2)){
  write_lines(out2[i], str_interp("${resdir_angsd_snapp.local}/hpd_heightfix2_tree${i}.nwk"))
}
# save summary df
write_csv(summary, str_interp("${resdir_angsd_snapp.local}/hpd_heightfix2_summary.csv"))

# treeanotation for 95%HPD trees
for(i in 1:length(out1)){
  system(str_interp("/home/aokidaisuke/miniconda3/envs/bioinfo/bin/treeannotator -heights median -burnin 50000 -target ${resdir_angsd_snapp.local}/hpd_heightfix2_tree${i}.nwk ${resdir_angsd_snapp.local}/snapp_heightfix2.trees ${resdir_angsd_snapp.local}/snapp_heightfix2_annotated${i}.trees"))
}

# treeannotation for maximum clade credibility tree
system(str_interp("/home/aokidaisuke/miniconda3/envs/bioinfo/bin/treeannotator -heights median -burnin 50000 ${resdir_angsd_snapp.local}/snapp_heightfix2.trees ${resdir_angsd_snapp.local}/snapp_heightfix2_annotatedMCC.trees"))

levels <- c("tomakomai", "kushiro", "russia", "aomori", "outgroup")

tree.df <- vector("list", length(out1))
for(i in 1:length(out2)){
  
  tree.df[[i]] <- read.beast(str_interp("${resdir_angsd_snapp.local}/snapp_heightfix2_annotated${i}.trees")) %>% as_tibble() 
  tree.df[[i]] <- 
    tree.df[[i]] %>% 
    mutate(
      percentage = summary$percentage[i],
      HPD_1 = map(height_0.95_HPD, .f = function(x)x[1]) %>% unlist %>% round(., 2),
      HPD_2 = map(height_0.95_HPD, .f = function(x)x[2]) %>% unlist %>% round(., 2),
      posterior.round = round(posterior, 2),
      height_median.round = signif(height_median, 3),
      posterior_95HPD = str_c(posterior.round, "\n", height_median.round, " [", HPD_1, ", ", HPD_2, "]")) %>% 
    mutate(label = if_else(str_detect(label, "Laterallus"), "outgroup", label) %>% factor(., levels = levels))
  tree.df[[i]] <- as.treedata(tree.df[[i]])
}

tree.mcc.df <- read.beast(str_interp("${resdir_angsd_snapp.local}/snapp_heightfix2_annotatedMCC.trees")) %>% as_tibble()

colv <- c(as.character(color("bright")(6))[c(2,4,1,3)], "#4D4D4D")
names(colv) <- levels

p.tree <- 
  ggdensitree(tree.df,aes(alpha = percentage*.01, color = percentage), size = 2) +
  geom_tippoint(aes(fill = label), pch = 21, size = 5)+
  geom_tiplab(aes(label=str_to_title(label)), size = 5, nudge_x = .3)+
  geom_nodelab(aes(x=-height_median, label=posterior_95HPD),hjust = 1, vjust = -.5, size=5) +
  geom_range(range = 'height_0.95_HPD', color = "grey66", alpha = .6, size = 2)+
  coord_cartesian(clip="off")+
  scale_color_viridis_c("Percentage", option = "B", direction = -1)+
  scale_fill_manual(values = colv)+
  guides(fill = "none", alpha = "none")+
  theme(
    plot.margin = unit(c(1, 3, 1, 1), "cm") ,
    legend.position = "bottom"
  )

p.tree
ggsave(str_interp("figures/revise/Fig2/snapp_hightfix2.png"))
ggsave(str_interp("figures/revise/Fig2/snapp_hightfix2.eps"))

