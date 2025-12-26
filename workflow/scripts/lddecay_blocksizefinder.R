#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript lddecay_blocksizefinder.R --ld LD.ld --acf input.acf.gz --treemix input.treemix --outdir outdir")
}

# Simple arg parser
get_arg <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 0 || (i+1) > length(args)) return(NA)
  args[i+1]
}

ld_path <- get_arg("--ld")
acf_path <- get_arg("--acf")
treemix_path <- get_arg("--treemix")
outdir <- get_arg("--outdir")

# Ensure output directory exists; write directly into it
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
plotdir <- outdir

# 1) Load ACF and compute retained site set excluding near-private alleles (as in original)
message("Reading ACF and building retained sites set...")

# Get header fields from acf
head <- system(
  paste("glactools view -P", shQuote(acf_path)),
  intern = TRUE
) %>% tail(1) %>% str_split("\t", simplify = TRUE) %>% as.character()

acf_df <- system(
  paste("glactools view", shQuote(acf_path)),
  intern = TRUE
) %>%
  str_split("\t", simplify = TRUE) %>%
  as_tibble(.name_repair = "minimal") %>%
  set_names(head) %>%
  select(!c(root, anc))

n.pop <- ncol(acf_df) - 3

# Mark zeros per population for ref and alt
ref <- acf_df %>%
  select(!c(`#chr`, coord, `REF,ALT`)) %>%
  mutate(across(everything(), ~str_extract(., "[0-9]+(?=,)"))) %>%
  set_names(paste0(names(.), "_ref")) %>%
  mutate(across(everything(), ~if_else(. == "0", "0", "1"))) %>%
  unite(comb_ref, everything(), sep = "", remove = TRUE) %>%
  mutate(countcomb_ref = str_count(comb_ref, "0"))

alt <- acf_df %>%
  select(!c(`#chr`, coord, `REF,ALT`)) %>%
  mutate(across(everything(), ~str_extract(., "(?<=,)[0-9]+"))) %>%
  set_names(paste0(names(.), "_alt")) %>%
  mutate(across(everything(), ~if_else(. == "0", "0", "1"))) %>%
  unite(comb_alt, everything(), sep = "", remove = TRUE) %>%
  mutate(countcomb_alt = str_count(comb_alt, "0"))

rm_df <- acf_df %>%
  bind_cols(ref, alt) %>%
  filter(countcomb_ref < (n.pop-1) & countcomb_alt < (n.pop-1))

rm_df_sub <- rm_df %>%
  select(!c(`#chr`, coord, `REF,ALT`, contains("comb"))) %>%
  mutate(across(everything(), ~str_remove(., ":[01]"))) %>%
  filter(if_all(everything(), ~ . != "0,0"))

input_treemix <- read_lines(treemix_path) %>%
  str_split(" ", simplify=TRUE) %>%
  as_tibble(.name_repair = "minimal") %>%
  { set_names(.[-1,], .[1,] %>% as.character()) }

check_match <- bind_cols(
  rm_df_sub %>% unite(comb_manual, everything(), sep ="/"),
  input_treemix %>% unite(comb_auto, everything(), sep ="/") %>% bind_rows(tibble(comb_auto=rep(NA_character_, nrow(rm_df_sub)-nrow(input_treemix))))
) %>% mutate(match = (comb_manual==comb_auto)) %>% filter(!match) %>% pull(match)
if(length(check_match)>0){
  stop("input.treemix and manually filtered sites do not match in the number of sites")
}

write_csv(rm_df %>% select(`#chr`, coord), file.path(outdir, "sites_input.treemix"))
sites <- read_csv(file.path(outdir, "sites_input.treemix"), show_col_types = FALSE) %>% unite(sites, everything(), sep=":") %>% pull(sites)

# 2) Load LD and filter to treemix sites
message("Reading LD and filtering to treemix sites...")
# ngsLD output (with --extend_out) has no header; the first row was being
# interpreted as column names, causing duplicate numeric names and missing
# variables like `site1`. Also, ngsLD may emit `-nan` which readr doesn't
# treat as NA by default. Read without header and rename the first 4 columns.
LD <- read_table(
    ld_path,
    col_names = FALSE,
    na = c("NaN", "nan", "-nan"),
    show_col_types = FALSE,
    col_types = cols(
      .default = col_skip(),
      X1 = col_character(),
      X2 = col_character(),
      X3 = col_double(),
      X4 = col_double()
    )
) %>%
  # Ensure standard names at least for the first 4 columns
  rename(site1 = X1, site2 = X2, dist = X3, r2 = X4) %>%
  # Keep only rows where both sites are in the retained set
  filter((site1 %in% sites) & (site2 %in% sites)) %>%
  # Drop rows with missing or non-finite values needed for fitting
  filter(is.finite(dist), is.finite(r2)) %>%
  arrange(dist)

if (nrow(LD) == 0) stop("No LD rows remain after site filtering. Check inputs.")

# 3) Smooth and compute half-life distance (distance where smoothed r2 drops to half its max)
message("Fitting LOESS to LD decay...")
lo <- loess(r2 ~ dist, data = LD, span = 0.2)
pred <- predict(lo, newdata = data.frame(dist = LD$dist))
ld_df <- tibble(distance = LD$dist, fpoints = pred)
half_val <- max(ld_df$fpoints, na.rm = TRUE) / 2
half_life <- ld_df %>% filter(fpoints > half_val) %>% summarise(max(distance)) %>% pull()
half_life <- ifelse(is.na(half_life), as.integer(median(LD$dist, na.rm = TRUE)), as.integer(round(half_life)))

write_csv(ld_df, file.path(plotdir, "ld.df.csv"))
write_csv(LD, file.path(plotdir, "ld_retained.csv"))
write_lines(half_life, file.path(plotdir, "block_size.txt"))

message(paste("Estimated block size =", half_life))

# 4) Plot
png(file.path(plotdir, "ld_plot.png"), height = 15, width = 15, units = "cm", res = 300)
ggplot() +
  geom_point(data = LD, aes(x=dist, y = r2), col = "#EE6677", alpha = 0.3)+
  geom_line(data = ld_df, aes(x=distance, y=fpoints))+
  xlab("Distance") + ylab(expression("LD "~(r^2)))+
  theme_bw()
invisible(dev.off())
