# SNAPP lognormal calibration helper ---------------------------------------
#
# This standalone snippet visualises the lognormal calibration prior that
# `snapp_prep.rb` injects into BEAST/SNAPP XML files. It is meant as a manual
# aid when you choose offset/mean/sd triplets for the constraint entries in
# the pipeline configuration; it is not sourced by any Snakemake rule.
#
# How to use:
#   1. Save any offset/mean/sd combination you want to inspect below.
#   2. Run the script in an interactive R session (`Rscript snapp_lognormal_helper.R`
#      or copy/paste into RStudio) to plot the density curve.
#   3. Adjust the numbers until the curve matches the calibration you intend to
#      encode (e.g., check that the offset matches the minimum age and the bulk
#      probability mass aligns with your fossil constraint).
#
# Construction notes:
#   * The numbers are interpreted exactly like `snapp_prep.rb` does: the first
#     entry is an offset (hard lower bound), the second is the real-space mean,
#     and the third the real-space standard deviation. BEAST uses a lognormal
#     parameterisation where `meanInRealSpace = true`, so we convert those
#     values to the log scale before calling `dlnorm`.
#   * The helper wraps the density in a function that clamps values below the
#     offset to zero, mirroring BEAST's behaviour for offset lognormals.
#
# Example parameters (replace with your own):
offset    <- 0.5   # Lower bound of the calibration (time units)
mean_real <- 11    # Mean of the distribution in real space
sd_real   <- 0.3   # Standard deviation in real space

# ---------------------------------------------------------------------------
# Convert the real-space mean/sd to mu/sigma on the log scale used by dlnorm.
sigma <- sqrt(log(1 + (sd_real^2 / mean_real^2)))
mu    <- log(mean_real) - 0.5 * sigma^2

# Offset lognormal density mirroring snapp_prep.rb / BEAST behaviour.
lognormal_density <- function(x) {
  y <- x - offset
  ifelse(y > 0, dlnorm(y, meanlog = mu, sdlog = sigma), 0)
}

# Choose a plotting window centred around the mean with extra tail coverage.
x_min <- offset
x_max <- mean_real + offset + 4 * sd_real

curve(
  lognormal_density,
  from = x_min,
  to = x_max,
  n = 500,
  xlab = "Time (same units as the calibration)",
  ylab = "Density",
  main = sprintf(
    "Lognormal calibration: offset = %.2f, mean = %.2f, sd = %.2f",
    offset, mean_real, sd_real
  )
)

# Visual guides: offset (minimum age) and the mean (offset + mean_real).
abline(v = offset, col = "firebrick", lty = 2)
abline(v = mean_real + offset, col = "steelblue", lty = 3)
legend(
  "topright",
  legend = c("Offset", "Mean"),
  col = c("firebrick", "steelblue"),
  lty = c(2, 3),
  bty = "n"
)

# Optional: print expected quantiles to the console for quick reference.
quantiles <- c(0.025, 0.5, 0.975)
cat("Quantiles (with offset applied):\n")
print(offset + qlnorm(quantiles, meanlog = mu, sdlog = sigma))
