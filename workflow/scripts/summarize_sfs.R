suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

# Helper to ensure directories exist
ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

# Extract structured metadata from file path
parse_components <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  n <- length(parts)
  list(
    site_filter = parts[n - 3],
    fold = parts[n - 2],
    pop = parts[n - 1]
  )
}

snakemake <- snakemake  # appease linters

sfs_files <- as.character(unlist(snakemake@input[["sfs"]]))
pest_files <- as.character(unlist(snakemake@input[["pest"]]))
contig_path <- snakemake@input[["contigs"]]

include_monomorphic <- isTRUE(as.logical(snakemake@params[["include_monomorphic"]]))
site_filters <- as.character(unlist(snakemake@params[["site_filters"]]))
fold_states <- as.character(unlist(snakemake@params[["fold_states"]]))
populations <- as.character(unlist(snakemake@params[["populations"]]))
min_contig_length <- as.numeric(snakemake@params[["min_contig_length"]])

if (length(sfs_files) == 0) {
  stop("No SFS files supplied to summarize_sfs.R")
}

# 1) SFS summaries ---------------------------------------------------------
read_sfs <- function(path) {
  comp <- parse_components(path)
  counts <- scan(path, what = numeric(), quiet = TRUE)
  tibble(
    file = path,
    site_filter = comp$site_filter,
    fold = comp$fold,
    pop = comp$pop,
    diff = seq_along(counts) - 1,
    frequency = counts
  )
}

sfs_df <- purrr::map_dfr(sfs_files, read_sfs) %>%
  mutate(
    site_filter = factor(site_filter, levels = unique(site_filters)),
    fold = factor(fold, levels = unique(fold_states)),
    pop = factor(pop, levels = populations)
  )

if (!include_monomorphic) {
  sfs_df <- sfs_df %>% filter(diff > 0)
}

# 2) Theta / Tajima summaries ---------------------------------------------
read_pest <- function(path) {
  comp <- parse_components(path)
  raw <- readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  colnames(raw) <- c(
    "idx_window", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL",
    "Tajima", "fuf", "fud", "fayh", "zeng", "nSites"
  )
  raw %>%
    mutate(
      across(
        c(
          "WinCenter", "tW", "tP", "tF", "tH",
          "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites"
        ),
        as.numeric
      )
    ) %>%
    mutate(
      file = path,
      site_filter = comp$site_filter,
      fold = comp$fold,
      pop = comp$pop
    )
}

pest_df <- purrr::map_dfr(pest_files, read_pest)

contigs <- readr::read_csv(contig_path, show_col_types = FALSE) %>%
  rename(Chr = accession, contig_length = length)

theta_df <- pest_df %>%
  left_join(contigs, by = "Chr") %>%
  mutate(
    contig_length = as.numeric(contig_length),
    site_filter = factor(site_filter, levels = unique(site_filters)),
    fold = factor(fold, levels = unique(fold_states)),
    pop = factor(pop, levels = populations)
  ) %>%
  filter(!is.na(pop))

theta_filtered <- theta_df %>%
  filter(!is.na(contig_length) & contig_length >= min_contig_length) %>%
  mutate(
    tW_per_site = if_else(nSites > 0, tW / nSites, NA_real_),
    tP_per_site = if_else(nSites > 0, tP / nSites, NA_real_)
  ) %>%
  filter(!is.na(tW_per_site), !is.na(Tajima), is.finite(tW_per_site), is.finite(Tajima))

# 3) Linear models ---------------------------------------------------------
fit_lm <- function(df, response, populations) {
  df_mod <- df %>%
    mutate(pop = factor(pop, levels = populations)) %>%
    filter(!is.na(pop), !is.na(contig_length), contig_length > 0)
  df_mod <- df_mod %>% filter(!is.na(.data[[response]]), is.finite(.data[[response]]))
  if (n_distinct(df_mod$pop) < 2) {
    return(NULL)
  }
  stats::lm(as.formula(paste(response, "~ pop + log(contig_length)")), data = df_mod)
}

model_tbl <- theta_filtered %>%
  group_by(site_filter, fold) %>%
  group_modify(~ {
    data <- .x
    tibble(
      metric = c("thetaW_per_site", "TajimaD"),
      response = c("tW_per_site", "Tajima"),
      fit = map(response, ~ fit_lm(data, .x, populations))
    )
  }) %>%
  ungroup() %>%
  filter(!map_lgl(fit, is.null)) %>%
  mutate(
    tidy = map(fit, ~ broom::tidy(.x, conf.int = TRUE)),
    glance = map(fit, broom::glance)
  )

model_summary <- model_tbl %>%
  select(site_filter, fold, metric, tidy) %>%
  unnest(tidy)

# 4) Plot ------------------------------------------------------------------
plot_data <- theta_filtered

plot_long <- plot_data %>%
  select(site_filter, fold, pop, tW_per_site, Tajima) %>%
  pivot_longer(
    cols = c(tW_per_site, Tajima),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      tW_per_site = "Theta W per site",
      Tajima = "Tajima's D"
    ),
    metric = factor(metric, levels = c("Theta W per site", "Tajima's D")),
    pop = factor(pop, levels = populations)
  ) %>%
  filter(!is.na(value), is.finite(value))

plot_obj <- ggplot(plot_long, aes(x = pop, y = value, fill = pop)) +
  geom_boxplot(width = 0.72, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.14, alpha = 0.45, size = 1.2, color = "black") +
  facet_grid(rows = vars(metric), cols = vars(site_filter, fold), scales = "free_y") +
  labs(
    x = "Population",
    y = "Window statistic",
    fill = "Population"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

n_site <- max(1, n_distinct(plot_data$site_filter))
n_fold <- max(1, n_distinct(plot_data$fold))
n_pop <- max(1, n_distinct(plot_long$pop))
plot_width <- max(8, 1.2 * n_pop * n_fold)
plot_height <- max(6, 3.2 * n_site * 2)

# 5) Write outputs ---------------------------------------------------------
ensure_parent(snakemake@output[["sfs_summary"]])
ensure_parent(snakemake@output[["theta_summary"]])
ensure_parent(snakemake@output[["lm_summary"]])
ensure_parent(snakemake@output[["plot"]])

readr::write_csv(sfs_df, snakemake@output[["sfs_summary"]])
readr::write_csv(theta_filtered, snakemake@output[["theta_summary"]])
if (nrow(model_summary) > 0) {
  readr::write_csv(model_summary, snakemake@output[["lm_summary"]])
} else {
  readr::write_csv(tibble(message = "Insufficient variation to fit models"),
                   snakemake@output[["lm_summary"]])
}

ggplot2::ggsave(
  filename = snakemake@output[["plot"]],
  plot = plot_obj,
  width = plot_width,
  height = plot_height,
  units = "in"
)
