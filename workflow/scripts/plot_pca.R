Sys.unsetenv("R_LIBS_USER")
Sys.unsetenv("R_PROFILE_USER")
Sys.unsetenv("R_ENVIRON_USER")

library(tidyverse)
library(magrittr)

safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)
palette_n <- function(n) {
  rep_len(safe_colorblind_palette, n)
}

# pcaplot function (from original script)
pcaplot <- function(xaxis, yaxis, .group, .pc, .e, .colpal, .labels) {
  plot_pc <- .pc %>%
    dplyr::select(all_of(c(xaxis, yaxis))) %>%
    mutate(group = .group) %>%
    set_colnames(c("xaxis", "yaxis", "group"))

  # Fix axis signs
  if (abs(min(plot_pc$xaxis)) > abs(max(plot_pc$xaxis))) plot_pc$xaxis <- -plot_pc$xaxis
  if (abs(min(plot_pc$yaxis)) > abs(max(plot_pc$yaxis))) plot_pc$yaxis <- -plot_pc$yaxis

  pct <- .e %>% filter(pc %in% c(xaxis, yaxis)) %>% pull(pct)
  labels <- str_to_upper(c(xaxis, yaxis))

  ellipse_df <- plot_pc %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    filter(n() >= 3) %>%
    ungroup()

  ggplot(plot_pc, aes(x = xaxis, y = yaxis, color = group)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point() +
    xlab(str_interp("${labels[1]} (${pct[1]}%)")) +
    ylab(str_interp("${labels[2]} (${pct[2]}%)")) +
    stat_ellipse(
      data = ellipse_df,
      level = 0.5,
      type = "norm",
      na.rm = TRUE
    ) +
    scale_color_manual(values = .colpal, labels = .labels) +
    theme_bw()
}

# Inputs
cov_mat <- as.matrix(read.table(snakemake@input[['cov']]))
eig_vals <- as.matrix(read.table(snakemake@input[['eig']]))

# Read bamlist (1-col TXT, no header)
bamlist <- read_delim(
  snakemake@input[['bamlist']],
  delim = "\t",
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  rename(path = 1) %>%
  mutate(sample = str_extract(path, "(?<=/)[^/]+(?=\\.bam)"))

# Samples metadata
samples <- read_tsv(snakemake@input[['samples']], show_col_types = FALSE)
populations <- snakemake@config[['populations']]
group_col <- snakemake@config[['group_col']]
axes <- snakemake@config[['pcangsd']][['axes_plot']]

# PCA data
pc <- eigen(cov_mat)$vectors %>%
  as.data.frame() %>%
  as_tibble() %>%
  set_colnames(paste0('pc', seq_len(ncol(.))))
e_df <- 
  eig_vals %>%
  as.matrix %>% 
  eigen %>% 
  .$values %>% 
  tibble(pc = str_c("pc",1:length(.)), values = .) %>% 
  mutate(pct = round(values/sum(values)*100, digits = 2))

# Merge metadata
df <- bamlist %>%
  left_join(samples, by = 'sample') %>%
  mutate(group_plot = factor(.data[[group_col]], levels = populations))

# Color palette
colpal <- as.character(palette_n(length(populations)))

# Generate plots for each PC combination
plots <- list(); idx <- 1
for (i in seq_along(axes)) {
  for (j in seq_along(axes)) {
    if (i < j) {
      plots[[idx]] <- pcaplot(
        axes[i], axes[j], df$group_plot,
        pc, e_df, colpal, populations
      )
      idx <- idx + 1
    }
  }
}
# saveRDS(list(plots = plots, df = df, pc = pc, e_df = e_df, populations = populations, axes = axes), "debug/plot_pca.RDS")

# Save all to single PDF
pdf(snakemake@output[[1]], width = 8, height = 6 * length(plots))
walk(plots, print)
dev.off()
