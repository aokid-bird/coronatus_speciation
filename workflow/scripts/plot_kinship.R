#!/usr/bin/env Rscript
Sys.unsetenv("R_LIBS_USER")
Sys.unsetenv("R_PROFILE_USER")
Sys.unsetenv("R_ENVIRON_USER")
library(tidyverse)
library(igraph)
library(glue)
library(readr)

safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)
palette_n <- function(n) {
  rep_len(safe_colorblind_palette, n)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: plot_kinship.R <kinfile> <datafile> <threshold> <popcol>")
}
kinfile <- args[1]
datafile <- args[2]
threshold <- as.numeric(args[3])
popcol <- args[4]

# extract output directory
resdir <- dirname(kinfile)

# read kinship matrix
kinmat <- 
  read_tsv(kinfile, show_col_types = FALSE) %>% 
  mutate(kinship = if_else(KING >= threshold, "related", "not-related"))

# plo1: kinship coefficient distribution
p <- ggplot(data = kinmat) +
  geom_point(aes(x = "a", y = KING, color = kinship)) + 
  geom_hline(yintercept = threshold) +
  theme_bw()
ggsave(glue("{resdir}/kin_KING_thr{threshold}.pdf"), plot = p)

# extract related pairs
kin.df <- kinmat %>% 
  filter(kinship == "related") %>% 
  select(ida, idb, KING) %>% 
  rename(from = ida, to = idb)

if (nrow(kin.df) == 0) {
  # create a csv
  write_csv(tibble(), glue("{resdir}/related_KING_thr{threshold}.csv"))
  pdf(file = glue("{resdir}/kinnet_KING_thr{threshold}.pdf"),
      width = 15, height = 15)
  plot.new()
  title(glue("No related pairs (KING >= {threshold})"))
  dev.off()
  quit(save="no")
}

# id lists
id <- unique(c(kin.df$from, kin.df$to))

# sample information
data.df <- read_tsv(datafile, show_col_types = FALSE)

kin.base.df <- data.df %>% 
  filter(sample %in% id) %>% 
  mutate(group = .[[popcol]], name = sample) %>%
  select(name, group)

# graph reconstruction and clustering
g <- graph_from_data_frame(kin.df, directed = FALSE, vertices = kin.base.df)
fgc <- cluster_fast_greedy(g)
grps <- groups(fgc) %>% lapply(as_tibble)

grps <- 
  lapply(seq_along(grps), function(x) mutate(grps[[x]], group = names(grps)[x])) %>% 
  bind_rows() %>% 
  rename(sample = value)

# join the kin results to base.df
tab.kin <- grps %>% rename(name = sample) %>% left_join(kin.base.df, by = "name")
write_csv(tab.kin, glue("{resdir}/related_KING_thr{threshold}.csv"))

# draw network
pdf(file = glue("{resdir}/kinnet_KING_thr{threshold}.pdf"),
    width = 15, height = 15, units = "cm", res = 300)

group_list <- sort(unique(V(g))$group)
group_colors <- setNames(palette_n(length(group_list)), group_list)
vertex_colors <- group_colors[V(g)$group]

plot(g,
     vertex.color = vertex_colors,
     edge.width = E(g)$KING * 20,
     edge.label = round(E(g)$KING, digits = 4),
     edge.label.color = "black",
     edge.label.cex = 2,
     vertex.size = 50,
     vertex.frame.width = 2,
     vertex.label.color = "black",
     vertex.label.cex = 1)
title(glue("Relatedness based on KING (threshold = {threshold})"))

legend("bottomright",
       legend = names(group_colors),
       col = group_colors,
       pch = 19,
       pt.cex = 2,
       cex = 0.9,
       bty = "n")

dev.off()
