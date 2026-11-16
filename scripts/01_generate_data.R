# =============================================================================
# scripts/01_generate_data.R
#
# Generate simulated survival data for the paper's two data-generation
# processes.
# Output per replication i:
#   data/sim_setting<k>/<scenario>/all/simdata_<i>/train.csv
#   data/sim_setting<k>/<scenario>/all/simdata_<i>/test.csv
#   data/sim_setting<k>/<scenario>/surv_th/simdata_<i>/times.csv
#   data/sim_setting<k>/<scenario>/surv_th/simdata_<i>/predth_test.csv
#   data/sim_setting<k>/<scenario>/surv_th/simdata_<i>/predth_alt_test.csv
#   data/sim_setting<k>/<scenario>/surv_th/simdata_<i>/predth_test_times.csv
#   data/sim_setting<k>/<scenario>/surv_th/simdata_<i>/predth_alt_test_times.csv
#   data/sim_setting<k>/<scenario>/surv_th/simdata_<i>/rmst_test_times.csv
#   data/sim_setting<k>/<scenario>/surv_th/simdata_<i>/rmst_alt_test_times.csv
# =============================================================================

library(here)
source(here("config.R"))

# =============================================================================
# USER SETTINGS 
# =============================================================================

SETTING <- 2    # 1 = Weibull/Biospear  |  2 = AFT/Friedman
S1_N    <- 2000 # paper uses 2000 or 20000
S2_N    <- 2000 # paper uses 2000 or 20000

# =============================================================================
# Setting 1 parameters 
# h(t|τ,X) = h0(t)·exp(α_τ·τ + Σ β_i1 f(X_i) + Σ β_i2 f(X_i)·τ)
# f(X) = (2X−1)^2,  α_τ = −0.288,  β_i1 = 0.288,  β_i2 = −0.9
# Weibull shape = 1 (exponential),  X_i ~ U[0,1]
# q = 20 biomarkers total; p = 10 active; censoring 20% or 50%
# =============================================================================
S1_PARAMS <- list(
  p     = 10,
  qmain = 10,
  qint  = 10,
  alpha = -0.288,
  beta  =  0.288,
  gamma = -0.9
)
S1_DATA_ROOT <- here("data", "sim_setting1")

# =============================================================================
# Setting 2 parameters 
# log T = m1(X) + τ·m2(X) + W,  W ~ Γ(2,1)
# Friedman random function generator, l=10 terms, subsets of size n_tot=5
# a2l ~ U[−0.2, 0.3],  X ~ N(0, I_20),  censoring ~20%
# =============================================================================
S2_PARAMS <- list(
  a_min   = -0.2,
  a_max   =  0.3,
  nvartot =  5,
  nvaract =  5,
  l       = 10,
  norm_ts = FALSE,
  cr      =  4.5
)
S2_DATA_ROOT <- here("data", "sim_setting2")

# =============================================================================
# Shared helpers — oracle survival at horizons + RMST
# =============================================================================

oracle_at_horizons <- function(Soracle, time_grid, horizons) {
  out <- matrix(1.0, nrow = nrow(Soracle), ncol = length(horizons))
  for (k in seq_along(horizons)) {
    idx <- which(time_grid <= horizons[[k]])
    if (length(idx) > 0) out[, k] <- as.numeric(Soracle[, max(idx)])
  }
  df <- as.data.frame(out)
  colnames(df) <- names(horizons)
  df
}

rmst_at_horizons <- function(Soracle, time_grid, horizons) {
  out <- matrix(0.0, nrow = nrow(Soracle), ncol = length(horizons))
  for (k in seq_along(horizons)) {
    tstar <- horizons[[k]]
    idx   <- which(time_grid <= tstar)
    if (length(idx) == 0) {
      out[, k] <- tstar
      next
    }
    i_lo    <- max(idx)
    i_hi    <- min(i_lo + 1L, length(time_grid))
    denom   <- time_grid[i_hi] - time_grid[i_lo]
    w       <- if (denom > 1e-10) (tstar - time_grid[i_lo]) / denom else 0
    s_tstar <- (1 - w) * Soracle[, i_lo] + w * Soracle[, i_hi]
    t_sub   <- c(0, time_grid[idx], tstar)
    S_sub   <- cbind(1, Soracle[, idx, drop = FALSE], s_tstar)
    dt      <- diff(t_sub)
    Savg    <- (S_sub[, -1, drop = FALSE] + S_sub[, -ncol(S_sub), drop = FALSE]) / 2
    out[, k] <- as.numeric(Savg %*% dt)
  }
  df <- as.data.frame(out)
  colnames(df) <- names(horizons)
  df
}

write_oracle_horizon_files <- function(surv_th_dir, times) {
  horizons <- c(P25 = as.numeric(times[["P25"]]),
                P5  = as.numeric(times[["P5"]]),
                P75 = as.numeric(times[["P75"]]))
  
  Soracle     <- as.matrix(read.csv(file.path(surv_th_dir, "predth_test.csv"),     row.names = 1))
  Soracle_alt <- as.matrix(read.csv(file.path(surv_th_dir, "predth_alt_test.csv"), row.names = 1))
  time_grid   <- as.numeric(sub("^X", "", colnames(Soracle)))
  
  write.csv(oracle_at_horizons(Soracle,     time_grid, horizons),
            file.path(surv_th_dir, "predth_test_times.csv"),     row.names = TRUE)
  write.csv(oracle_at_horizons(Soracle_alt, time_grid, horizons),
            file.path(surv_th_dir, "predth_alt_test_times.csv"), row.names = TRUE)
  write.csv(rmst_at_horizons(Soracle,       time_grid, horizons),
            file.path(surv_th_dir, "rmst_test_times.csv"),       row.names = TRUE)
  write.csv(rmst_at_horizons(Soracle_alt,   time_grid, horizons),
            file.path(surv_th_dir, "rmst_alt_test_times.csv"),   row.names = TRUE)
}

# =============================================================================
# Run data generation
# =============================================================================
if (SETTING == 1) {
  
  source(here("src", "data", "simulate_setting1.R"))
  dir.create(S1_DATA_ROOT, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("Setting 1 | n=%d | %d replications", S1_N, length(ITERS)))
  
  for (i in ITERS) {
    message(sprintf("  [Setting 1] replication %d / %d", i, max(ITERS)))
    simulate_setting1(
      i         = i,
      n         = S1_N,
      p         = S1_PARAMS$p,
      qmain     = S1_PARAMS$qmain,
      qint      = S1_PARAMS$qint,
      alpha     = S1_PARAMS$alpha,
      beta      = S1_PARAMS$beta,
      gamma     = S1_PARAMS$gamma,
      data_root = paste0(S1_DATA_ROOT, "/")
    )
    scenario_dir <- list.dirs(S1_DATA_ROOT, recursive = FALSE)[1]
    surv_th_dir  <- file.path(scenario_dir, "surv_th", paste0("simdata_", i))
    times        <- read.csv(file.path(surv_th_dir, "times.csv"))
    write_oracle_horizon_files(surv_th_dir, times)
  }
  message("Setting 1 data generation complete.")
  
} else if (SETTING == 2) {
  
  source(here("src", "data", "simulate_setting2.R"))
  dir.create(S2_DATA_ROOT, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("Setting 2 | n=%d | %d replications", S2_N, length(ITERS)))
  
  for (i in ITERS) {
    message(sprintf("  [Setting 2] replication %d / %d", i, max(ITERS)))
    simulate_setting2(
      num_data  = i,
      n         = S2_N,
      a_min     = S2_PARAMS$a_min,
      a_max     = S2_PARAMS$a_max,
      nvartot   = S2_PARAMS$nvartot,
      nvaract   = S2_PARAMS$nvaract,
      l         = S2_PARAMS$l,
      norm_ts   = S2_PARAMS$norm_ts,
      cr        = S2_PARAMS$cr,
      data_root = paste0(S2_DATA_ROOT, "/")
    )
    scenario_dir <- list.dirs(S2_DATA_ROOT, recursive = FALSE)[1]
    surv_th_dir  <- file.path(scenario_dir, "surv_th", paste0("simdata_", i))
    times        <- read.csv(file.path(surv_th_dir, "times.csv"))
    write_oracle_horizon_files(surv_th_dir, times)
  }
  message("Setting 2 data generation complete.")
  
} else {
  stop("SETTING must be 1 or 2.")
}

