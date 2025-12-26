#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript treemix_optm_eval.R --runs_dir DIR --max_m N --out_prefix PREFIX [--reps N] [--plotting_funcs RPATH] [--optm]")
}

get_arg <- function(flag, default=NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || (i+1) > length(args)) return(default)
  args[i+1]
}

runs_dir <- get_arg("--runs_dir")
max_m <- as.integer(get_arg("--max_m"))
out_prefix <- get_arg("--out_prefix")
reps_arg <- suppressWarnings(as.integer(get_arg("--reps", NA)))
plotting_funcs <- get_arg("--plotting_funcs", NA)
use_optm <- any(args == "--optm")

if (!is.na(plotting_funcs) && nzchar(plotting_funcs) && file.exists(plotting_funcs)) {
  try(source(plotting_funcs), silent = TRUE)
} else if (!is.na(plotting_funcs) && nzchar(plotting_funcs)) {
  message(sprintf("Plotting funcs not found at %s; skipping variance explained.", plotting_funcs))
}

# Gather .llik files
paths <- list.files(runs_dir, pattern = "\\.llik$", full.names = TRUE, recursive = TRUE)
if (length(paths) == 0) stop("No .llik files found for evaluation")

# Replicate original extraction: for each file, read lines and extract
# the numeric value following 'events: ' from each line; then coerce into
# a 2-column matrix by row (V1, V2) and bind with file metadata.
vals <- paths %>%
  lapply(read_lines) %>%
  lapply(stringr::str_extract, "(?<=events: )[0-9.-]+(?= )") %>%
  unlist(use.names = FALSE) %>%
  suppressWarnings(as.numeric(.))

if (length(vals) %% 2 != 0) {
  warning("Unexpected number of extracted values from .llik files; results may be misaligned.")
}

df_base <- matrix(vals, ncol = 2, byrow = TRUE) %>%
  as_tibble(.name_repair = "minimal") %>%
  `colnames<-`(c("V1", "V2"))

llik_df <- df_base %>%
  mutate(path = paths) %>%
  mutate(edge = suppressWarnings(as.numeric(stringr::str_extract(path, "(?<=_e)[0-9]+(?=_)"))),
         rep  = suppressWarnings(as.numeric(stringr::str_extract(path, "(?<=_o)[0-9]+")))) %>%
  mutate(V1 = as.numeric(V1), V2 = as.numeric(V2))

# Derive reps vector: prefer --reps; else infer from files
reps_vec <- if (!is.na(reps_arg) && reps_arg > 0) seq_len(reps_arg) else sort(unique(llik_df$rep))

#===========================#
# 5-1) Variance Explained  #
#===========================#
var_df <- NULL
if (exists("calcVarExplain")) {
  # Loop over edges (0..max_m) and reps, mirroring original script
  var_accum <- list()
  for (e in 0:max_m) {
    for (o in reps_vec) {
      prefix <- file.path(runs_dir, sprintf("treemix_e%d_o%d", e, o))
      vfile <- paste0(prefix, ".vertices.gz")
      efile <- paste0(prefix, ".edges.gz")
      if (file.exists(vfile) && file.exists(efile) && file.info(vfile)$size > 0 && file.info(efile)$size > 0) {
        ve <- try(calcVarExplain(prefix), silent = TRUE)
        if (!inherits(ve, "try-error")) {
          var_accum[[length(var_accum)+1]] <- tibble(edge=e, o=o, VarExplain=as.numeric(ve[[2]]))
        }
      }
    }
  }
  if (length(var_accum) > 0) {
    var_df <- bind_rows(var_accum)
    write_csv(var_df, paste0(out_prefix, "_variance.csv"))
  }
} else {
  message("calcVarExplain() not found (plotting_funcs not sourced?); skipping variance explained.")
}

#===========================#
# 5-2) Combine with lnL     #
#===========================#
params_df <- llik_df %>%
  # Join variance explained by edge/rep (o), as in original
  { if (!is.null(var_df)) left_join(., var_df, by=c("edge" = "edge", "rep" = "o")) else mutate(., VarExplain = NA_real_) } %>%
  mutate(diff = ifelse(is.finite(V1) & is.finite(V2), V2 - V1, NA_real_))

write_csv(params_df, paste0(out_prefix, "_parameters.csv"))
write_csv(params_df, paste0(out_prefix, "_runs.csv"))

# Summary per edge: lnL mean/sd/se and VarExplain mean/sd/se
summary_df <- params_df %>%
  group_by(edge) %>%
  summarise(
    mean    = mean(V2, na.rm = TRUE),
    sd      = sd(V2, na.rm = TRUE),
    se      = sd / sqrt(sum(!is.na(V2))),
    mean.var = mean(VarExplain, na.rm = TRUE),
    sd.var   = sd(VarExplain, na.rm = TRUE),
    se.var   = sd.var / sqrt(sum(!is.na(VarExplain))),
    .groups = 'drop'
  )

write_csv(summary_df, paste0(out_prefix, "_summary.csv"))

#===========================#
# OptM (Evanno-style)      #
#===========================#
if (use_optm) {
  suppressWarnings(suppressMessages({
    try(library(OptM), silent = TRUE)
  }))
  if ("package:OptM" %in% search() || requireNamespace("OptM", quietly = TRUE)) {
    try({
      optm <- OptM::optM(runs_dir, tsv = paste0(out_prefix, "_optm.tsv"))
      # Create PDF externally and request plots without internal file writing
      pdf_file <- paste0(out_prefix, "_optm.pdf")
      grDevices::pdf(pdf_file)
      plots <- try(OptM::plot_optM(optm, plot = FALSE), silent = TRUE)
      if (!inherits(plots, "try-error")) {
        # Handle common return types: ggplot object, list of plots, or grob
        if (inherits(plots, "gg") || inherits(plots, "ggplot")) {
          print(plots)
        } else if (is.list(plots)) {
          for (p in plots) {
            if (inherits(p, "gg") || inherits(p, "ggplot")) {
              print(p)
            } else if (inherits(p, "grob") || inherits(p, "gTree") || inherits(p, "gList")) {
              grid::grid.newpage(); grid::grid.draw(p)
            } else {
              # Fallback: try to print whatever it is
              try(print(p), silent = TRUE)
            }
          }
        } else if (inherits(plots, "grob") || inherits(plots, "gTree") || inherits(plots, "gList")) {
          grid::grid.newpage(); grid::grid.draw(plots)
        } else {
          # Fallback to printing
          try(print(plots), silent = TRUE)
        }
      }
      grDevices::dev.off()
    }, silent = TRUE)
  } else {
    message("OptM not available in this environment; skipping.")
  }
}
