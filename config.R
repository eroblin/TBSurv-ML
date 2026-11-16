# =============================================================================
# config.R
# =============================================================================

# Number of simulation sets
ITERS <- 1:10

# Column labels for the three time horizons (25th, 50th, 75th percentile of
# the marginal event-time distribution, computed per simulation replicate)
TIME_LABELS <- c("P25", "P5", "P75")

# Number of bootstrap resamples for RATE/AUTOC confidence intervals
RATE_BOOTSTRAP_R <- 200
