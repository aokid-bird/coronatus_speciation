#==================================================#
#### Figure 2 and related supplementary figures ####
#==================================================#
# This script intends to plot figures related to phylogenetic analyses
# RAxML runs were conducted with different subsets of samples, specifically, those retaining or
# excluding Kanto wintering samples and those retaining or excluding invariant sites.
# Figure 2 will be presented only with representative ones, while the others will be saved as 
# supplementary materials.

#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
library(treeio)
library(ggtree)
library(khroma)
library(colorblindr)
library(ggnewscale)
source("code/dataprep.R")
source("../functions/ngs_functions.R")
source("../functions/cluster_functions.R")
source("~/treemix-1.13/src/plotting_funcs.R")
source("../functions/general_functions.R")
source("code/ggplot_treemix.R")

#==============================#
#### Directories & base dfs ####
#==============================#
##### Directories #####
resdir_angsd_raxml.local1 <- "res/angsd_raxml_rlrmv_mindepthind2_minind0.5"
resdir_angsd_raxml.local2 <- "res/angsd_raxml_rlrmv_mindepthind2_minind0.8"
resdir_angsd_treemix.local <- "res/angsd_treemix_mindepthind2"
rundir_ngsdist.local <- "res/ngsdist_mindepthind2"
resdir_bwa.local <- "res/bwa"
resdir_bwa_sra.local <- "res/bwa_sra"
figdir <- "figures/revise"
# clusters
resdir_bwa_sra.clust <- str_interp("${dir_clust}/res/bwa_sra")
resdir_bwa.clust <- str_interp("${dir_clust}/${resdir_bwa.local}")

##### base dfs #####
subdir <- 
  expand.grid(category = c("all", "breed"),
              outgrp   = c("_wout", ""),
              variant = c("", "_onlyvar")) %>% 
  unite("subdir", category:variant, sep="") %>% 
  pull(subdir)

summary.df <- 
  read_csv("res/bcf/summary_df.csv") %>% 
  filter(call.depth == 4) %>% 
  mutate(X1 = str_c(dir_clust, "/", bamlist.fullpath))

#==================#
#### 2) Treemix ####
#==================#
subdir <- 
  expand.grid(category = c("breed"),
              outgrp   = c("", "_wout"),
              pop = c("pop", "indiv")) %>% 
  mutate(subdir = str_c(category, outgrp)) %>% 
  select(subdir, pop)

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


ggtreemix <- resid <- vector("list", length = 6)
ggtreemix <- rep(list(ggtreemix), 4)
resid <- rep(list(resid), 4)

for(k in 1:nrow(subdir)){
  
  tmp.subdir <- subdir$subdir[k]
  tmp.pop    <- subdir$pop[k]
  
  resdir_treemix <- str_interp("${resdir_angsd_treemix.local}/${tmp.subdir}/treemix_${tmp.pop}")
  
  summary_trmx.df <- 
    read_csv(str_interp("${resdir_treemix}/result_summary.csv"))
  
  llk_trmx.df <-
    read_csv(str_interp("${resdir_treemix}/parameters.csv"))
  
  optm <- read_tsv(str_interp("${resdir_treemix}/optm.tsv"))
  
  best <- 
    llk_trmx.df %>% 
    group_by(edge) %>% 
    mutate(max_score = max(V2)) %>% 
    ungroup() %>% 
    filter(V2==max_score) %>% 
    group_by(edge) %>% 
    slice(1) %>% 
    ungroup() 
  
  ##### 2-1) Evaluation plots #####
  # plotting evaluation results
  # log likelihood
  
  llk.p_multi <- 
    ggplot() + 
    geom_point(data = llk_trmx.df, aes(x = edge, y = V2)) + 
    geom_point(data = summary_trmx.df, aes(x = edge, y = mean), 
               color = as.character(color("high contrast")(3)[3]))+
    geom_line(data = summary_trmx.df, aes(x = edge, y = mean), 
              color = as.character(color("high contrast")(3)[3]))+
    xlab("Number of migration edge") + ylab("Log likelihood") + 
    ggtitle("a) Mean Log-likelihood") + 
    theme_bw()+
    theme(legend.position = "non",
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)) 
  
  varexp.p_multi <- 
    ggplot() + 
    geom_point(data = llk_trmx.df, aes(x = edge, y = VarExplain))+
    geom_point(data = summary_trmx.df, aes(x = edge, y = mean.var), 
               color = as.character(color("high contrast")(3)[3])) + 
    geom_line(data = summary_trmx.df, aes(x = edge, y = mean.var), 
              color = as.character(color("high contrast")(3)[3]))  + 
    geom_hline(yintercept = 0.998, lty = 2)+
    xlab("Number of migration edge") + ylab("Variance explained")+ 
    ggtitle("b) Variance explained") +
    theme_bw()+
    scale_y_continuous(n.breaks = 10)+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)) 
  
  ##### 2-2) Residual plots ##### 
  # resid <- vector("list", 6)
  
  for(i in 1:6){
    bestrep <- best$rep[i]
    edge <- i-1
    
    resid[[k]][[i]]<-plot_resid(stem = str_interp("${resdir_treemix}/treemix_e${edge}_o${bestrep}"), 
                           pop_order = str_interp("${resdir_treemix}/popordcol"),
                           cex = 0.4, usemax = F)
    
  }
  
  rng <- lapply(resid[[k]], range) %>% lapply(matrix, ncol = 2) %>% lapply(as_tibble) %>% bind_rows 
  min.rng <- min(rng$V1) - abs(min(rng$V1)*0.02)
  max.rng <- max(rng$V2) + abs(max(rng$V1)*0.02)
  
  crep <- best %>% filter(edge == 1) %>% pull(rep)
  plot_tree(stem = str_interp("${resdir_treemix}/treemix_e1_o${crep}"))
  validation
  png(str_interp("figures/revise/treemix_sup/treemix_validation_${tmp.subdir}_${tmp.pop}.png"))
     OptM::plot_optM(optm)
  dev.off()
  postscript(str_interp("figures/revise/treemix_sup/treemix_validation_${tmp.subdir}_${tmp.pop}.eps"),
             width = 10, height = 20)
    OptM::plot_optM(optm)
  dev.off()
  
  llk.p_multi / varexp.p_multi
  ggsave(str_interp("figures/revise/treemix_sup/treemix_validation2_${tmp.subdir}_${tmp.pop}.png"))
  ggsave(str_interp("figures/revise/treemix_sup/treemix_validation2_${tmp.subdir}_${tmp.pop}.eps"),
         width = 10, height = 20, units = "cm")
  
  # residual plots with multiple edges
  alpha <- c("a)", "b)", "c)", "d)", "e)", "f)")
  
  png(str_interp("figures/revise/treemix_sup/treemix_resid_${tmp.subdir}_${tmp.pop}.png"),
      width = 20, height = 20, units = "cm", res = 300)
  par(mfrow = c(3, 2), mar = c(1,1,4,1), oma = c(2,2,2,3))
  for(i in 1:6){
    bestrep <- best$rep[i]
    edge <- i-1
    lab <- alpha[i]
    VE <- signif(best$VarExplain[i], 3)
    LLK <- signif(best$V2[i],3)
    resid[[i]]<-plot_resid(stem = str_interp("${resdir_treemix}/treemix_e${edge}_o${bestrep}"), 
                           pop_order = str_interp("${resdir_treemix}/popordcol"),
                           cex = 0.4, usemax = F, max = max.rng, min = min.rng)
    title(str_interp("${lab} m = ${edge} (LogLike = ${LLK}, VarExp = ${VE})"))
  }
  dev.off()
  postscript(str_interp("figures/revise/treemix_sup/treemix_resid_${tmp.subdir}_${tmp.pop}.eps"),
                        width = 20, height = 20)
  par(mfrow = c(3, 2), mar = c(1,1,4,1), oma = c(2,2,2,3))
  for(i in 1:6){
    bestrep <- best$rep[i]
    edge <- i-1
    lab <- alpha[i]
    VE <- signif(best$VarExplain[i], 3)
    LLK <- signif(best$V2[i],3)
    resid[[i]]<-plot_resid(stem = str_interp("${resdir_treemix}/treemix_e${edge}_o${bestrep}"),
                           pop_order = str_interp("${resdir_treemix}/popordcol"),
                           cex = 0.4, usemax = F, max = max.rng, min = min.rng)
    title(str_interp("${lab} m = ${edge} (LogLike = ${LLK}, VarExp = ${VE})"))
  }
  dev.off()
  
  tree plots with mutiple edges
  par(mfrow = c(3, 2), mar = c(1,1,4,1), oma = c(2,2,2,3))
  for(i in 1:6){
    
    bestrep <- best$rep[i]
    edge <- i-1
    lab <- alpha[i]
    VE <- signif(best$VarExplain[i], 3)
    LLK <- signif(best$V2[i],3)
    obj <- read_treemix(stem = str_interp("${resdir_treemix}/treemix_e${edge}_o${bestrep}"))
    levels <- c("tomakomai", "kushiro", "russia", "aomori", "kanto", "Laja")
    
    if(tmp.pop == "pop"){
      
      pop.df <-
        bamlist.df %>% 
        mutate(group.label = factor(group.label, levels = levels)) %>% 
        select(group.label) %>% 
        distinct() %>% 
        mutate(sample = group.label)
      
    }else if(tmp.pop == "indiv"){
      
      pop.df <-
        bamlist.df %>%
        mutate(group.label = factor(group.label, levels = levels)) %>% 
        select(sample, group.label) 
    }
    
    pop.df %<>% mutate(group.label)
    
    # color
    color.palette <- 
      tibble(group_reg = levels,
             group_color = c(as.character(color("bright")(6))[c(2,4,1,3,5)], "#4D4D4D")) %>% 
      filter(group_reg %in% unique(pop.df$group.label))
    
    ggtreemix[[k]][[i]] <- plot_treemix(obj=obj,
                                   pop.df=pop.df,
                                   .col.pal = color.palette$group_color) 
    ggtreemix[[k]][[i]] <- 
      ggtreemix[[k]][[i]] + 
      ggtitle(str_interp("Migration edge = ${i-1}"))+
      theme(
        plot.margin = unit(c(1, 1, 1, 1), "cm")  
      )+
      coord_cartesian(clip = "off") + 
        theme_treemix()
    
    # plot_tree(stem = str_interp("${resdir_treemix}/treemix_e${edge}_o${bestrep}"))
    # title(str_interp("${lab} m = ${edge} (LogLike = ${LLK}, VarExp = ${VE})"))
  }
  
  ggtreemix[[k]][[1]]+ggtreemix[[k]][[2]]+ggtreemix[[k]][[3]]+
    ggtreemix[[k]][[4]]+ggtreemix[[k]][[5]]+ggtreemix[[k]][[6]]+
    plot_layout(nrow = 2)
  # ggsave(str_interp("figures/revise/Fig2/treemix_all_${tmp.subdir}_${tmp.pop}.png"))
}

ggtreemix[[2]][[2]]
ggsave("figures/revise/Fig2/treemix_pop.eps",device=cairo_ps);ggsave("figures/revise/Fig2/treemix_pop.png")
ggtreemix[[4]][[2]]
ggsave("figures/revise/Fig2/treemix_indiv.eps", device=cairo_ps);ggsave("figures/revise/Fig2/treemix_indiv.png")

#===========================================#
##### 5-3) TreeMix - multiple best runs #####

figdir_treemix <- "figures/revise/treemix_sup"
k = 2

ggtreemix[[k]][[1]]+ggtreemix[[k]][[2]]+ggtreemix[[k]][[3]]+
  ggtreemix[[k]][[4]]+ggtreemix[[k]][[5]]+ggtreemix[[k]][[6]]+
  plot_layout(nrow = 3)+
  plot_annotation(tag_levels = "a", tag_suffix = ")")

ggsave(str_interp("${figdir_treemix}/best_tree_rep.png"), width = 20, height = 30, units = "cm", dpi = 500)
ggsave(str_interp("${figdir_treemix}/best_tree_rep.eps"), device = cairo_ps, width = 20, height = 30, units = "cm")

#===========================================#
##### 5-3) TreeMix - multiple best runs #####
alpha <- c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)")

resdir_treemix <- "res/angsd_treemix_mindepthind2/breed_wout/treemix_pop"

postscript(str_interp("${figdir_treemix}/mig_1_rep.eps"), width = 20, height = 35)
par(mfrow = c(5, 2), mar = c(1,1,4,1), oma = c(2,2,2,3))
for(i in 1:10){
  llk_trmx.df <-
    read_csv(str_interp("${resdir_treemix}/parameters.csv"))
  crep <- i
  .edge <- 1
  lab <- alpha[i]
  VE <- signif(llk_trmx.df %>% filter(rep == crep & edge == .edge) %>% pull(VarExplain), 3)
  LLK <- signif(llk_trmx.df %>% filter(rep == crep & edge == .edge) %>% pull(V2),3)
  plot_tree(stem = str_interp("${resdir_treemix}/treemix_e${edge}_o${crep}"))
  title(str_interp("${lab} Replicate ${crep} (LogLike = ${LLK}, VarExp = ${VE})"))
}
dev.off()
