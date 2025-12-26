#=========================#
#### Diversity Indices ####
#=========================#
# This script intends to compare and plot group-level diversity indices, including tajima's D and thetas
# from the results of angsd_sfs.R.
# Comparisons are conducted using linear modelling.

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
library(lmerTest)
library(lme4)
library(DHARMa)
library(broom)
library(ggsignif)
library(igraph)
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
resdir_angsd_sfs.local <- "res/angsd_sfs_mindepthind2"
resdir_bwa_sra.local <- "res/bwa_sra"
resdir_divind.local <- "res/diversity_indices_mindepthind2"

mkdir(resdir_divind.local)

##### 1-2) Cluster #####
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")
resdir_angsd_sfs.clust <- str_interp("${dir_clust}/${resdir_angsd_sfs.local}")
resdir_bwa_sra.clust <- str_interp("${dir_clust}/res/bwa_sra")

#=========================#
#### 2) dataframe prep ####
#=========================#
bamlist.df <-
  data.df %>%
  filter(!remove) %>%
  filter(group_reg != "kanto_w") %>%
  mutate(ref = "data.df") %>%
  select(sample, group_2, ref) %>%
  rename(group = group_2)

subdir <-
  bamlist.df %>%
  pull(group) %>%
  unique() %>%
  as.character()

ls.sfs <- list.files(resdir_angsd_sfs.local, recursive = TRUE, pattern = "*fold.sfs", full.names = TRUE)
ls.sum <- list.files(resdir_angsd_sfs.local, recursive = TRUE, pattern = "*thetas.idx.pestPG", full.names = TRUE)
pat.dir <- str_c("?<=", resdir_angsd_sfs.local, "/")

# SFSs
sfs_glob.df <-
  ls.sfs %>%
  lapply(., read_delim, delim = " ", col_name = FALSE) %>%
  lapply(., function(x){select(x, !colnames(x)[ncol(x)])}) %>%
  lapply(., function(x){
    t(x) %>%
      as_tibble %>%
      mutate(diff = row_number()-1) %>%
      rename(frequency = V1) %>%
      select(diff, frequency)}) %>%
  mapply(function(x,y){
    mutate(x, path = y)
  }, x = ., y = ls.sfs, SIMPLIFY = FALSE) %>%
  bind_rows %>%
  mutate(pop = str_extract(path, str_interp("(${pat.dir}).*(?=/gl_)"))%>% factor(., levels = c("insular", "continental")),
         fold = str_extract(path, "(?<=gl_intersect_).*(?=.sfs)") %>% str_extract("fold|unfold"),
         link = str_extract(path, "unlinked") %>% if_else(is.na(.), "linked", .))

# diversity indices summary
summary.df <-
  ls.sum %>%
  lapply(., read_delim) %>%
  lapply(., magrittr::set_colnames, c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)",
                                      "Chr",	"WinCenter",	"tW",	"tP",	"tF",	"tH",	"tL",
                                      "Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")) %>%
  mapply(function(x,y){
    mutate(x, path = y)
  }, x = ., y = ls.sum, SIMPLIFY = FALSE) %>%
  bind_rows %>%
  mutate(pop = str_extract(path, str_interp("(${pat.dir}).*(?=/gl_)"))%>% factor(., levels = c("insular", "continental")),
         fold = str_extract(path, "(?<=gl_intersect_).*(?=.thetas)") %>% str_extract("fold|unfold"),
         link = str_extract(path, "unlinked") %>% if_else(is.na(.), "linked", .))

write_csv(sfs_glob.df, str_interp("${resdir_divind.local}/sfs_df.csv"))
write_csv(summary.df, str_interp("${resdir_divind.local}/divind_summary_df.csv"))

#===========================#
#### 3) Statistical test ####
#===========================#
#===================#
##### 3-1) Prep #####
# contig length data
contig <- read_csv("/mnt/y/Studies/data/WGS/Laterallus_jamaicensis_coturniculus/contig_length.csv")

# filter out contigs smaller than 2M base
contig.list <- contig %>% filter(length > 2000000) %>% pull(accession)

theta.df <-
  summary.df %>%
  filter(Chr %in% contig.list) %>%
  left_join(., contig, by = c("Chr" = "accession")) %>%
  mutate(tW_persite = tW/nSites,
         tP_persite = tP/nSites) %>%
  filter(link == "linked") %>%
  group_by(fold) %>%
  nest()

#========================#
##### 3-2) Modelling #####
# Theta
fun <- function(d, base, resp.var){
  l.pop <- as.character(unique(d$pop))
  sub.pop <- l.pop[l.pop != base]
  d.mod <-
    d %>%
    mutate(pop = factor(pop, levels = c(base, sub.pop))) %>%
    select(c(all_of(resp.var), pop, length)) %>%
    magrittr::set_colnames(c("resp", "pop", "length"))
  d.mod %<>% filter(!is.infinite(resp))
  fit <- lm(resp ~ pop + log(length),
            data = d.mod)
  return(fit)
}
#### ----fun end ####
res.fit <-
  theta.df %>%
  mutate(fit_insular_thetaW   = map(data, fun, base = "insular", resp.var = "tW_persite"),
         # fit_thetaPi  = map(data, fun, base = "insular", resp.var = "tP_persite"),
         fit_insular_tajimaD  = map(data, fun, base = "insular", resp.var = "Tajima"),
         # fit_fufs     = map(data, fun, base = "insular", resp.var = "fuf"),
  ) %>%
  pivot_longer(cols = starts_with("fit_"),
               names_to = c("base1", "index"),
               values_to = "fit",
               names_pattern = "fit_(.*)_(.*)") %>%
  mutate(est = map(fit, tidy, conf.int = TRUE)) %>%
  unnest(est) %>%
  mutate(base2 = if_else(!(term %in% c("(Intercept)", "log(length)")), str_remove(term, "pop"), NA)) %>%
  # model evaluation by DHARMa
  mutate(simresid = map(fit, simulateResiduals, plot = F))
# save RDS
saveRDS(res.fit, str_interp("${resdir_divind.local}/res.fit.rds"))
res.fit <- readRDS(str_interp("${resdir_divind.local}/res.fit.rds"))

# save csv
res.fit %>% filter(fold == "fold") %>%
  select(fold, base1, index, term, estimate:base2) %>%
  write_csv(str_interp("${resdir_divind.local}/res.fit.csv"))


#===================#
#### 4) plotting ####
#===================#
##### 4-1) SFS #####
ggplot(data=sfs_glob.df, aes(x = diff, y = log(frequency), fill = pop, group = fold)) +
  geom_col() +
  facet_grid(cols = vars(fold, pop), rows = vars(link), scales = "free_y") +
  ggtitle("SFS including 0")
# ggsave("figures/sfs/SFS1D_with0.png")

sfs_glob.df %>%
  filter(diff != 0) %>%
  ggplot(data=, aes(x = diff, y = frequency, fill = pop, group = fold)) +
  geom_col() +
  facet_grid(cols = vars(fold, pop), rows = vars(link), scales = "free") +
  ggtitle("SFS excluding 0")
