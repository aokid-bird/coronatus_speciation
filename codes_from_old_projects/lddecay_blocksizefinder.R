#================#
#### LD decay ####
#================#
# This script intends to fit LD-decay curve to the LD calcualted by ngsLD and to identify
# linkage block size for TreeMix

#============#
#### Args ####
#============#
args <- commandArgs(trailingOnly = TRUE)
resdir_angsd_global.local <-args[1]
tmp.subdir <- args[2]
resdir_treemix <- args[3]

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
library(plyr)
source("../functions/ngs_functions.R")
source("../functions/cluster_functions.R")
source("../functions/general_functions.R")
source("code/dataprep.R")

rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange

# loaded directories
#### local
# datadir_local <- "../data/NGS/migseq_200906"
# refpath_local <- "../data/WGS/Laterallus_jamaicensis_coturniculus/JAKCOX01.fasta"
# refdir_local <-  "../data/WGS/Laterallus_jamaicensis_coturniculus"
#### cluster
# dir_clust <- str_interp("/lfs/user/${currdir}")
# datadir_clust <- str_replace(datadir_local, "..", "/lfs/user")
# refdir_clust  <- str_replace(refdir_local, "..", "/lfs/user")
# refpath_clust <- str_replace(refpath_local, "..", "/lfs/user")
# dirclust <- "/lfs/user/${currdir}"
data.df
ref.df

#======================#
#### 1) Directories ####
#======================#
# resdir_angsd_global.local <- "res/angsd_global"
# tmp.subdir <- "breed"
# resdir_treemix <- "res/angsd_treemix/breed_wout/treemix_indiv"
# if no lddecay directory created yet
mkdir(str_interp("${resdir_treemix}/lddecay"))

#====================#
#### 2) Data Prep ####
#====================#
#==================================#
##### 2-1) genotype locus data #####
geno <-
  system(str_interp("zcat ${resdir_angsd_global.local}/${tmp.subdir}/gl_global.beagle.gz"), intern = TRUE) %>%
  str_split("\t", simplify = TRUE) %>%
  as.data.frame

geno.colname <-
  geno[1,]

geno %<>% slice(-1)

geno %<>% select(V1:V3)

#================================#
##### 2-2) treemix site data #####
head <-
  str_interp("glactools view -P ${resdir_treemix}/input.acf.gz") %>%
  system(intern = TRUE) %>%
  .[length(.)] %>%
  str_split("\t", simplify = TRUE) %>%
  as.character

acf.df <-
  str_interp("glactools view ${resdir_treemix}/input.acf.gz") %>%
  system(intern = TRUE) %>%
  str_split("\t", simplify = TRUE) %>%
  as_tibble %>%
  set_colnames(head) %>%
  select(!c(root,anc))

n.pop <- ncol(acf.df) - 3 #

# count any ref/alt sites that contains 0 for more than n-1 populations out of n
# which means private alleles
ref <- 
  acf.df %>%
  select(!c(`#chr`, coord, `REF,ALT`)) %>%
  mutate_all(~str_extract(., "[0-9]+(?=,)")) %>%
  set_colnames(str_c(colnames(.), "_ref")) %>%
  mutate_all(~if_else(.=="0","0","1")) %>%
  unite(comb_ref, everything(), sep = "", remove = TRUE) %>%
  mutate(countcomb_ref = str_count(comb_ref, "0"))

alt <- acf.df %>%
  select(!c(`#chr`, coord, `REF,ALT`)) %>%
  mutate_all(~str_extract(., "(?<=,)[0-9]+")) %>%
  set_colnames(str_c(colnames(.), "_alt")) %>%
  mutate_all(~if_else(.=="0","0","1")) %>%
  unite(comb_alt, everything(), sep = "", remove = TRUE) %>%
  mutate(countcomb_alt = str_count(comb_alt, "0"))

# remove those sites with private alleles
rm.df <-
  acf.df %>%
  bind_cols(., ref, alt) %>%
  filter(countcomb_ref < (n.pop-1) & countcomb_alt < (n.pop-1))

# select the populations to be checked and also
# filter out sites that have "0,0"
rm.df_sub <-
  select(rm.df, !c(`#chr`, coord, `REF,ALT`,contains("comb"))) %>%
  mutate_all(~str_remove(., ":[01]")) %>%
  filter_all(all_vars(. != "0,0"))

# load the treemix input to be comparable with the rm.df_sub
input.treemix <- str_interp("${resdir_treemix}/input.treemix")
tree.input <- 
  read_lines(input.treemix) %>%
  str_split(" ", simplify=TRUE) %>%
  as_tibble %>%
  set_colnames(.[1,]) %>% slice(-1)

# check all the sites manually retained are the same as treemix input file
check.match <- 
  bind_cols(
  rm.df_sub %>% unite(comb_manual, everything(), sep ="/"),
  tree.input %>% unite(comb_auto, everything(), sep ="/") %>% bind_rows(tibble(comb_auto=rep(NA_character_, nrow(rm.df_sub)-nrow(tree.input))))
) %>%
  mutate(match = (comb_manual==comb_auto)) %>%
  filter(!match) %>% 
  pull(match)
if(length(check.match)>0){
  stop("input.treemix and manually filtered sites do not match in the number of sites")
}
# TRUE
nrow(rm.df)

# Output the retained site names
rm.df %>%
  select(`#chr`, coord) %>%
  write_csv(str_interp("${resdir_treemix}/sites_input.treemix"))
# reload the site files
sites <-
  read_csv(str_interp("${resdir_treemix}/sites_input.treemix")) %>%
  unite(sites, everything(), sep=":") %>%
  pull(sites)

# LD results with retaining only the sites that TreeMix used
LD <-
  read_table(str_interp("${resdir_angsd_global.local}/${tmp.subdir}/LD.ld")) %>%
  filter((site1 %in% sites) & (site2 %in% sites)) %>% 
  arrange(dist)

# parameters
Cstart <- c(C=0.1)
n=length(unique(as.character(geno.colname)))-3
modelC <- nls

#=====================#
#### 2) LD fitting ####
#=====================#
modelC <- LDdecay(df = LD, Cstart=Cstart, fit = TRUE)

new.rho <- summary(modelC)$parameters[1]

fpoints <- LDdecay(df = LD, rho = new.rho, fit = FALSE)

ld.df <- tibble(distance = LD$dist, .fpoints = fpoints)

half_life <- ld.df %>% filter(.fpoints > max(.fpoints)/2) %>% pull(distance) %>% max

# Half life = 22

write_csv(ld.df, str_interp("${resdir_treemix}/lddecay/ld.df.csv"))
write_csv(LD, str_interp("${resdir_treemix}/lddecay/ld_retained.csv"))
write_lines(half_life, str_interp("${resdir_treemix}/lddecay/block_size.txt"))

#===================#
#### 3) Plotting ####
#===================#
y.label <- max(ld.df$.fpoints)/2
x.label <- half_life + (500/10)
p.ld <- 
  ggplot() +
  geom_point(data = LD, aes(x=dist, y = r2), col = "#EE6677", alpha = 0.3)+
  geom_line(data = ld.df, aes(x=distance, y=fpoints))+
  xlim(c(0,500))+
  geom_vline(xintercept = half_life, color = "#004488", lty = 2)+
  geom_hline(yintercept = y.label, color = "#004488", lty = 2)+
  geom_label(aes(x=x.label, y=y.label), label = str_interp("block size = ${half_life}"))+
  xlab("Distance")+
  ylab(expression("LD "~(r^2)))+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

png(str_interp("${resdir_treemix}/lddecay/ld_plot.png"), height = 15, width = 15, units = "cm", res = 300)
print(p.ld)
dev.off()