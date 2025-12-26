.libPaths(c("/home/aokid/R/x86_64-redhat-linux-gnu-library/4.1", .libPaths()))

library(tidyverse)
library(magrittr)
library(khroma)

# Inputs
qopt_files <- snakemake@input[['qopt']]
log_files  <- snakemake@input[['logs']]
beagle_f   <- snakemake@input[['beagle']]
bamlist_f  <- snakemake@input[['bamlist']]
samples    <- read_tsv(snakemake@input[['samples']], show_col_types = FALSE)
populations <- snakemake@config[['populations']]
group_col  <- snakemake@config[['group_col']]
maxK       <- snakemake@config[['ngsadmix']][['maxK']]

# Read log-likelihoods
df_like <- map_dfr(log_files, function(f) {
  lines <- readLines(f)
  like_line <- grep('^b', lines, value = TRUE)
  tibble(
    K      = as.integer(str_extract(basename(f), '(?<=K)\\d+')),
    run    = str_remove(basename(f), '\\..+$'),
    like   = as.numeric(str_extract(like_line, '(?<=like=)[0-9.-]+'))
  )
})

# DeltaK plot
df_delta <- df_like %>%
  group_by(K) %>%
  summarise(deltaK = mean(abs(like)) / sd(abs(like))) %>%
  filter(K > 1)

p1 <- ggplot(df_delta, aes(x = K, y = deltaK)) +
  geom_line() +
  theme_bw() +
  labs(x = 'K', y = 'delta K')

ggsave(snakemake@output[['delta']], p1, width = 6, height = 4)

# Select best run per K
best <- df_like %>% group_by(K) %>% filter(like == max(like)) %>% slice(1)

# Derive sample order from Beagle header (3 cols + 3 per individual)
# Read the first line directly to avoid name-repair side effects
hdr_line <- readLines(beagle_f, n = 1)
tokens <- strsplit(hdr_line, "\t", fixed = TRUE)[[1]]
if (length(tokens) < 4) {
  stop('Beagle file seems malformed: fewer than 4 columns in header: ', beagle_f)
}
triplet_names <- tokens[4:length(tokens)]
# ANGSD beagle repeats each sample name 3 times; take first of each triplet
if (length(triplet_names) %% 3 == 0) {
  idx <- seq(1, length(triplet_names), by = 3)
  sample_order <- triplet_names[idx]
} else {
  warning('Beagle header after first 3 columns is not a multiple of 3; inferring sample names by unique()')
  sample_order <- unique(triplet_names)
}

# Derive sample IDs from bamlist (paths per line)
bam_paths <- readr::read_lines(bamlist_f)
bam_paths <- bam_paths[nzchar(bam_paths)]
if (length(bam_paths) == 0) {
  stop('Bamlist appears empty: ', bamlist_f)
}
bam_samples <- basename(bam_paths) %>% str_remove('\\.(bam|cram)$')
if (length(bam_samples) != length(sample_order)) {
  stop('Bamlist entries (', length(bam_samples), ') do not match number of individuals inferred from Beagle header (',
       length(sample_order), '). Check for mismatches in bamlist and ANGSD inputs.')
}
plot_sample_levels <- bam_samples

# Ensure samples metadata includes at least 'sample' and add Beagle IDs for traceability
if (!('sample' %in% colnames(samples))) {
  stop("samples TSV must include a 'sample' column: ", snakemake@input[['samples']])
}
samples_joined <- tibble(
    sample = bam_samples,
    beagle_id = sample_order
  ) %>%
  left_join(samples, by = 'sample')

# Load Qopt and merge with ordered samples
df_list <- list()
for (i in seq_len(nrow(best))) {
  k  <- best$K[i]
  run <- best$run[i]
  f  <- file.path(dirname(qopt_files[1]), paste0(run, '.qopt'))
  tmp <- read.table(f) %>% as_tibble()
  colnames(tmp) <- paste0('V', seq_len(ncol(tmp)))
  if (nrow(tmp) != nrow(samples_joined)) {
    stop('Q-matrix rows (', nrow(tmp), ') do not match Beagle-derived sample count (',
         nrow(samples_joined), ') for ', basename(f), '. This usually means the metadata sample list and the Beagle/NGSadmix inputs differ. Check inclusion of outgroups and sample IDs.')
  }
  df_list[[i]] <- bind_cols(samples_joined, tmp) %>% mutate(K = k)
}
df_admix <- bind_rows(df_list)

# Pivot and plot
long <- df_admix %>%
  pivot_longer(cols = starts_with('V'), names_to = 'popGroup', values_to = 'prob') %>%
  mutate(
    Population = factor(.data[[group_col]], levels = populations),
    K = factor(K, levels = 1:maxK)
  ) %>%
  filter(K != "1")

pal <- as.character(color('muted')(maxK))

# Use bamlist-derived sample labels while preserving Beagle order
p2 <- ggplot(long, aes(x = factor(sample, levels = plot_sample_levels), y = prob, fill = popGroup)) +
  geom_col(color = 'grey', size = 0.1) +
  facet_grid(rows = vars(K), cols = vars(Population), scales = 'free', space = 'free') +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = 'sample')

ggsave(snakemake@output[['plots']], p2, width = 8, height = 4 * maxK)
