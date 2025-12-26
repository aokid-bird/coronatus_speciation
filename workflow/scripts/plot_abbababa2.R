suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
})

khroma_loaded <- requireNamespace("khroma", quietly = TRUE)
if (!khroma_loaded) {
  khroma_lib <- snakemake@params[["khroma_lib"]]
  if (!is.null(khroma_lib) && nzchar(khroma_lib)) {
    .libPaths(c(khroma_lib, .libPaths()))
    khroma_loaded <- requireNamespace("khroma", quietly = TRUE)
  }
}
if (!khroma_loaded) {
  stop("Package 'khroma' not available. Install it or set abbababa2.khroma_lib in config.")
}

khroma_colour <- function(name) {
  get("colour", asNamespace("khroma"))(name)
}

df <- readr::read_csv(snakemake@input[["summary"]], show_col_types = FALSE)

pop_codes <- snakemake@config[["populations"]]
pop_labels <- snakemake@config[["population_labels"]]
map_label <- function(code) {
  if (!is.null(pop_labels) && code %in% names(pop_labels)) {
    lbl <- pop_labels[[code]]
    if (!is.null(lbl) && !is.na(lbl) && nzchar(lbl)) {
      return(lbl)
    }
  }
  code
}

population_levels <- vapply(pop_codes, map_label, character(1))

plot_df <- df %>%
  filter(file == "Uncorrected") %>%
  mutate(
    H1 = factor(H1, levels = population_levels),
    H2 = factor(H2, levels = population_levels),
    H3 = factor(H3, levels = population_levels)
  ) %>%
  arrange(H1, H2, H3) %>%
  mutate(signif = label)

pal <- as.character(khroma_colour("bright")(2))

p <- ggplot(plot_df, aes(x = pair, y = D)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_linerange(aes(ymin = D - 3 * sd, ymax = D + 3 * sd), color = pal[1], linewidth = 0.8) +
  geom_linerange(aes(ymin = D - sd, ymax = D + sd), color = pal[2], linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(aes(label = signif), nudge_x = 0.15, nudge_y = 0.007) +
  coord_flip() +
  theme_bw() +
  labs(x = "Population pair", y = "D-statistic") +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11)
  )

out_dir <- dirname(snakemake@output[["pdf"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(snakemake@output[["pdf"]], plot = p, width = 8, height = 8)
