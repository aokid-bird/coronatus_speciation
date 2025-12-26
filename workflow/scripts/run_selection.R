# scripts/run_selection.R
args <- commandArgs(trailingOnly = TRUE)
library(R.utils)
library(tidyverse)
library(RcppCNPy)

selection_file <- args[1]
sites_file <- args[2]
out_csv <- args[3]
out_png <- args[4]

site <- readLines(sites_file)
D <- npyLoad(selection_file)
p <- pchisq(D, 1, lower.tail = FALSE)

qqchi <- function(x, ...) {
  lambda <- round(median(x) / qchisq(0.5, 1), 2)
  qqplot(qchisq((1:length(x) - 0.5) / (length(x)), 1), x, ylab = "Observed", xlab = "Expected", ...)
  abline(0, 1, col = 2, lwd = 2)
  abline(0, 1)
  legend("topleft", paste("lambda=", lambda))
}

png(out_png, width = 1500, height = 1500)
qqchi(D)
dev.off()

selection <- tibble(site = as.integer(site), p = p[, 1]) %>%
  mutate(selection = if_else(p < (0.05 / length(p)), 1, 0))
write_csv(selection, out_csv)
