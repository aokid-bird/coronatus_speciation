#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(ggnewscale)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript treemix_make_plots.R --runs_dir DIR --mode pop|indiv --samples TSV --group_col COL --bamlist PATH --max_m N --reps N --outdir DIR --plotting_funcs RPATH --ggplot_treemix RPATH")
}

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || (i + 1) > length(args)) {
    return(default)
  }
  args[i + 1]
}

runs_dir <- get_arg("--runs_dir")
mode <- get_arg("--mode", "pop")
samples_tsv <- get_arg("--samples")
group_col <- get_arg("--group_col")
bamlist_path <- get_arg("--bamlist")
max_m <- suppressWarnings(as.integer(get_arg("--max_m", "6")))
reps <- suppressWarnings(as.integer(get_arg("--reps", "10")))
outdir <- get_arg("--outdir")
plotting_funcs <- get_arg("--plotting_funcs")
ggplot_treemix_path <- get_arg("--ggplot_treemix")

if (!dir.exists(runs_dir)) stop("runs_dir not found: ", runs_dir)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load khroma only from a specific lib if needed, without altering global .libPaths
khroma_loaded <- FALSE
try({
  suppressPackageStartupMessages(library(khroma, quietly = TRUE, warn.conflicts = FALSE))
  khroma_loaded <- TRUE
}, silent = TRUE)
if (!khroma_loaded) {
  khroma_lib <- Sys.getenv("KHROMA_LIB")
  if (nzchar(khroma_lib)) {
    try({
      suppressPackageStartupMessages(library(khroma, lib.loc = khroma_lib, quietly = TRUE, warn.conflicts = FALSE))
      khroma_loaded <- TRUE
    }, silent = TRUE)
  }
}
if (!requireNamespace("khroma", quietly = TRUE)) {
  stop("Package 'khroma' not found. Install it in the conda env or set KHROMA_LIB to its library path.")
}
khroma_colour <- function(name) {
  get("colour", asNamespace("khroma"))(name)
}

# Source helper scripts
if (!is.null(plotting_funcs) && nzchar(plotting_funcs) && file.exists(plotting_funcs)) {
  source(plotting_funcs)
}
if (!is.null(ggplot_treemix_path) && nzchar(ggplot_treemix_path) && file.exists(ggplot_treemix_path)) {
  source(ggplot_treemix_path)
}

# Load eval outputs
eval_runs_csv <- file.path(runs_dir, "eval_runs.csv")
eval_summary_csv <- file.path(runs_dir, "eval_summary.csv")
if (!file.exists(eval_runs_csv)) {
  # fallback to prefixed file
  csvs <- list.files(runs_dir, pattern = "_runs\\.csv$", full.names = TRUE)
  if (length(csvs) >= 1) eval_runs_csv <- csvs[[1]]
}
if (!file.exists(eval_summary_csv)) {
  csvs <- list.files(runs_dir, pattern = "_summary\\.csv$", full.names = TRUE)
  if (length(csvs) >= 1) eval_summary_csv <- csvs[[1]]
}

llk_trmx <- readr::read_csv(eval_runs_csv, show_col_types = FALSE)
summary_trmx <- readr::read_csv(eval_summary_csv, show_col_types = FALSE)

# Determine best replicate per edge by maximum log-likelihood (V2)
best <- llk_trmx %>%
  group_by(edge) %>%
  mutate(max_score = max(V2, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(V2 == max_score) %>%
  group_by(edge) %>%
  slice(1) %>%
  ungroup()

# Helper to open a PDF device safely on clusters
save_plot_pdf <- function(filename, width = 10, height = 8) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(filename, width = width, height = height, paper = "special")
}

# 1) Validation plots: lnL and Var. explained
llk_p <- ggplot() +
  geom_point(data = llk_trmx, aes(x = edge, y = V2)) +
  geom_point(data = summary_trmx, aes(x = edge, y = mean), color = as.character(khroma_colour("high contrast")(3)[3])) +
  geom_line(data = summary_trmx, aes(x = edge, y = mean), color = as.character(khroma_colour("high contrast")(3)[3])) +
  xlab("Number of migration edge") +
  ylab("Log likelihood") +
  ggtitle("a) Mean Log-likelihood") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10))

var_p <- ggplot() +
  geom_point(data = llk_trmx, aes(x = edge, y = VarExplain)) +
  geom_point(data = summary_trmx, aes(x = edge, y = mean.var), color = as.character(khroma_colour("high contrast")(3)[3])) +
  geom_line(data = summary_trmx, aes(x = edge, y = mean.var), color = as.character(khroma_colour("high contrast")(3)[3])) +
  geom_hline(yintercept = 0.998, lty = 2) +
  xlab("Number of migration edge") +
  ylab("Variance explained") +
  ggtitle("b) Variance explained") +
  theme_bw() +
  scale_y_continuous(n.breaks = 10) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))

val_pdf <- file.path(outdir, "treemix_validation2.pdf")
save_plot_pdf(val_pdf, width = 10 / 2.54, height = 20 / 2.54)
print(llk_p / var_p)
grDevices::dev.off()

# 2) Residuals across best runs for m=0..max_m
# Create pop_order file by reading one run's cov header to get order
pick_edge <- if (0 %in% best$edge) 0 else min(best$edge, na.rm = TRUE)
pick_rep <- best %>%
  filter(edge == pick_edge) %>%
  slice(1) %>%
  pull(rep)
stem0 <- file.path(runs_dir, sprintf("treemix_e%d_o%d", pick_edge, pick_rep))
cov_path <- paste0(stem0, ".cov.gz")
if (!file.exists(cov_path)) stop("cov.gz not found for run: ", stem0)
cov_df <- suppressWarnings(read.table(gzfile(cov_path), as.is = TRUE, head = TRUE, quote = "", comment.char = ""))
pop_names <- names(cov_df)

# Build pop.df based on mode and samples metadata for coloring in ggplot trees
samples_df <- readr::read_tsv(samples_tsv, show_col_types = FALSE)
if (!("sample" %in% colnames(samples_df))) stop("samples TSV must include a 'sample' column")
if (!(group_col %in% colnames(samples_df))) stop(paste0("samples TSV must include the group column: ", group_col))
group_map <- samples_df %>%
  select(sample, !!sym(group_col)) %>%
  rename(group.label = !!sym(group_col))

if (mode == "pop") {
  pop_df <- tibble(sample = pop_names, group.label = pop_names)
} else {
  pop_df <- tibble(sample = pop_names) %>% left_join(group_map, by = "sample")
  # Fill NAs (e.g., outgroups/taxa labels) with their own name
  pop_df <- pop_df %>% mutate(group.label = if_else(is.na(group.label), sample, group.label))
}

# khroma 'bright' palette in original order
uniq_groups <- unique(pop_df$group.label)
pal <- khroma_colour("bright")(max(3, length(uniq_groups)))
group_colors <- setNames(pal[seq_along(uniq_groups)], uniq_groups)

# Write popordcol file used by treemix_plotting_funcs.R
popord_file <- file.path(runs_dir, "popordcol")
readr::write_delim(tibble(V1 = pop_names, V2 = group_colors[pop_df$group.label] %>% unname()), popord_file, delim = " ", col_names = FALSE)

# Residuals range across best runs for consistent color scaling
resid_list <- vector("list", max_m + 1)
for (i in 0:max_m) {
  bestrep <- best %>%
    filter(edge == i) %>%
    slice(1) %>%
    pull(rep)
  if (length(bestrep) == 0 || is.na(bestrep)) next
  stem <- file.path(runs_dir, sprintf("treemix_e%d_o%d", i, bestrep))
  # get residual matrix via plot_resid (returns matrix invisibly)
  resid_list[[i + 1]] <- try(plot_resid(stem = stem, pop_order = popord_file, cex = 0.4, usemax = FALSE), silent = TRUE)
}
rng_vals <- resid_list %>%
  purrr::keep(~ !inherits(., "try-error")) %>%
  lapply(range) %>%
  lapply(matrix, ncol = 2) %>%
  lapply(as_tibble) %>%
  bind_rows()
min_rng <- min(rng_vals$V1, na.rm = TRUE) - abs(min(rng_vals$V1, na.rm = TRUE) * 0.02)
max_rng <- max(rng_vals$V2, na.rm = TRUE) + abs(max(rng_vals$V2, na.rm = TRUE) * 0.02)

resid_pdf <- file.path(outdir, sprintf("treemix_resid_%s.pdf", mode))
save_plot_pdf(resid_pdf, width = 20 / 2.54, height = 20 / 2.54)
par(mfrow = c(3, 2), mar = c(1, 1, 4, 1), oma = c(2, 2, 2, 3))
for (i in 0:min(5, max_m)) {
  bestrep <- best %>%
    filter(edge == i) %>%
    slice(1) %>%
    pull(rep)
  if (length(bestrep) == 0 || is.na(bestrep)) next
  stem <- file.path(runs_dir, sprintf("treemix_e%d_o%d", i, bestrep))
  ve <- best %>%
    filter(edge == i) %>%
    pull(VarExplain) %>%
    signif(3)
  llk <- best %>%
    filter(edge == i) %>%
    pull(V2) %>%
    signif(3)
  plot_resid(stem = stem, pop_order = popord_file, cex = 0.4, usemax = FALSE, max = max_rng, min = min_rng)
  title(sprintf("%s) m = %d (LogLike = %s, VarExp = %s)", letters[i + 1], i, llk, ve))
}
grDevices::dev.off()

# 3) Tree plots for multiple edges (0..5 default grid)

trees_pdf <- file.path(outdir, sprintf("treemix_trees_%s.pdf", mode))

plots <- list()
for (i in 0:min(5, max_m)) {
  bestrep <- best %>%
    filter(edge == i) %>%
    slice(1) %>%
    pull(rep)
  if (length(bestrep) == 0 || is.na(bestrep)) next
  stem <- file.path(runs_dir, sprintf("treemix_e%d_o%d", i, bestrep))
  obj <- read_treemix(stem = stem)
  pd <- pop_df %>% select(sample, group.label)
  col_vec <- group_colors[pd$group.label] %>% unname()
  plt <- plot_treemix(obj = obj, pop.df = pd, .col.pal = col_vec) +
    ggtitle(sprintf("Migration edge = %d", i)) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    coord_cartesian(clip = "off") +
    theme_treemix()
  plots[[length(plots) + 1]] <- plt
}

if (length(plots) > 0) {
  combo <- wrap_plots(plots, nrow = 2)
  save_plot_pdf(trees_pdf, width = 24 / 2.54, height = 16 / 2.54)
  print(combo)
  grDevices::dev.off()
}

message("TreeMix plotting completed. Outputs in ", outdir)
