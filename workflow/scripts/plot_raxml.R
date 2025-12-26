.libPaths(c("/home/aokid/R/x86_64-redhat-linux-gnu-library/4.1", .libPaths()))

library(tidyverse)
library(magrittr)
library(patchwork)
library(khroma)
library(ggnewscale)
library(treeio)
library(ggtree)

tmp.bamlist <- read_delim(snakemake@input[['bamlist']], delim = " ", col_names = FALSE)
data.df <- read_tsv(snakemake@input[["samples"]])
geno <- snakemake@input[["geno"]]
raxsup <- snakemake@input[["raxsup"]]
group_col <- snakemake@params[["group_col"]]
populations <- snakemake@params[["populations"]]
# info data for each sample
levels <- populations
df.info <-
    tmp.bamlist %>% 
    mutate(sample = str_remove(X1, ".*/") %>% str_remove(".bam") %>% str_remove("_slice.*")) %>% 
    left_join(., data.df, by = "sample") %>% 
    mutate(pop = .[[group_col]]) %>%
    select(sample, pop) %>%
    mutate(pop = factor(pop, levels = levels))

# missing data information for the used sites
geno <- 
    system(str_interp("zcat ${geno}"), intern = TRUE) %>% 
    str_split("\t", simplify = TRUE) %>% 
    as_tibble %>% 
    rename(scaf = V1, pos = V2) %>% 
    unite("scaf_pos", c("scaf", "pos"), sep = "_")
# remove empty columns
geno %<>%  
    select(!colnames(geno)[ncol(geno)])

df.info <- 
    geno %>% 
    t %>% 
    as_tibble %>% 
    set_colnames(.[1,]) %>% 
    slice(-1) %>% 
    unite("concat", everything(), sep = "") %>% 
    bind_cols(df.info, .) %>% 
    mutate(countN = str_count(concat, "NN"),
            missingprop = countN/(str_width(concat)/2)) %>% 
    select(!concat)

# color palette
color.palette <- 
    tibble(pop = levels,
            group_color = c(as.character(color("bright")(length(populations))))) %>% 
    filter(pop %in% unique(df.info$pop))

#### 4-2) Tree plotting ####
# load the tree result with bootstrap
tree <- read.tree(raxsup)

# add information dataset and create a base plot
tmp <- 
    ggtree(tree) %<+% 
    df.info + 
    geom_tippoint(aes(fill=pop), size = 4, shape = 21, color = "black") + 
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color = "#1A1A1A", size = 2, shape = 16) +
    # geom_tiplab(aes(color=group_reg), size = 2.5, shape = 21,hjust = -0.4) + 
    scale_fill_manual(breaks = color.palette$pop,
                        values = color.palette$group_color, name = "Population") +
    scale_color_manual(breaks = color.palette$pop,
                        values = color.palette$group_color, name = "Population") +
    geom_treescale()+
    theme(
        plot.margin = unit(c(1, 3, 1, 1), "cm"),
        legend.position = "bottom"
    )+
    coord_cartesian(clip = "off")

p.raxml <- 
    tmp + 
    geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
                family = "helvetica", vjust=-.5, hjust = 1.2)+
    ggtitle("RAxML-ng")

#### 4-3) Tree plotting with QC info ####
# Tree with missing proportion + mean coverage data

df.sub <- 
    df.info %>% 
    select(sample, missingprop) %>% 
    rename(prop.miss = missingprop) %>% 
    mutate(sample = factor(sample, levels = tree$tip.label)) %>% 
    arrange(sample) %>% 
    mutate(sample=as.character(sample))

#saveRDS(list(df.sub = df.sub, tree = tree), "debug/tree.RDS")

p.raxml_wQC <-
    tmp + new_scale_fill()+
    geom_facet(panel = "Mean Missing Proportion", data = df.sub, geom = geom_col,
                mapping = aes(x = prop.miss, fill = prop.miss), orientation = "y") +
    theme_tree2()+
    scale_fill_viridis_c(option="D")

ggsave(snakemake@output[["plots"]], p.raxml_wQC, width = 8, height = 8)
