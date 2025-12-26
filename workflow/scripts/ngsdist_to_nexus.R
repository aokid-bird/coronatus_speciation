#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3){
  stop("Usage: ngsdist_to_nexus.R <ngsdist_matrix> <out_nexus> <labels_yes_no>")
}

mat_path <- args[1]
out_path <- args[2]
labels_flag <- args[3] # "yes" or "no"

source("workflow/scripts/ngs_functions.R")

# ngsDist outputs start with two header lines
mat <- readr::read_tsv(mat_path, skip = 2, col_names = FALSE) %>%
  dplyr::rename(ind = X1)

dist2nexus(df=mat, dryrun = FALSE, dir = out_path, labels = labels_flag, triangle = "both")

cat(sprintf("Wrote NEXUS: %s\n", out_path))

