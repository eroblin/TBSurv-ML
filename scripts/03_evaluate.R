# =============================================================================
# scripts/03_evaluate.R
#
# Compute all evaluation measures for one model on simulated data.
# Set SETTING, MODEL, and MATCH_TYPE below, then source from the repo root.
#
# Output:
#   results/sim_setting<k>/<MODEL>/resAll_CI_match<MATCH_TYPE>.csv
# =============================================================================

library(here)
source(here("config.R"))
source(here("src", "metrics", "run_scenario.R"))

# =============================================================================
# USER SETTINGS 
# =============================================================================
SETTING    <- 2         # 1 = Weibull/Biospear  |  2 = AFT/Friedman
MODEL      <- "ALASSO"   # "ALASSO" | "InteractionForest" | "CSF" | "CoxCC" | "CoxTime"
MATCH_TYPE <- "predben"  # "predben" | "var"

# =============================================================================
# Paths and file names
# =============================================================================
if (SETTING == 1) {
  data_root  <- list.dirs(here("data", "sim_setting1"), recursive = FALSE)[1]
  res_root   <- here("results", "sim_setting1")
} else if (SETTING == 2) {
  data_root  <- list.dirs(here("data", "sim_setting2"), recursive = FALSE)[1]
  res_root   <- here("results", "sim_setting2")
} else {
  stop("SETTING must be 1 or 2.")
}

sim_details <- basename(data_root)
pred_root   <- file.path(res_root, MODEL, sim_details)
oracle_root <- file.path(data_root, "surv_th")
output_dir  <- file.path(res_root, MODEL)

pred_basename            <- if (MODEL == "CSF") "pred_horizons_sf"  else "pred_horizons"
pred_alt_basename        <- "pred_horizons_alt"
rmst_basename            <- if (MODEL == "CSF") "rmst_horizons_sf"  else "rmst_horizons"
rmst_alt_basename        <- "rmst_horizons_alt"
oracle_basename          <- "predth_test_times"
oracle_alt_basename      <- "predth_alt_test_times"
rmst_oracle_basename     <- "rmst_test_times"
rmst_oracle_alt_basename <- "rmst_alt_test_times"

# =============================================================================
# Run evaluation
# =============================================================================
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf("Evaluating | setting %d | model: %s | match: %s | scenario: %s",
                SETTING, MODEL, MATCH_TYPE, sim_details))

run_scenario(
  scenario_name            = sprintf("Setting%d", SETTING),
  data_root                = data_root,
  pred_root                = pred_root,
  oracle_root              = oracle_root,
  iters                    = ITERS,
  match_type               = MATCH_TYPE,
  pred_basename            = pred_basename,
  pred_alt_basename        = pred_alt_basename,
  rmst_basename            = rmst_basename,
  rmst_alt_basename        = rmst_alt_basename,
  oracle_basename          = oracle_basename,
  oracle_alt_basename      = oracle_alt_basename,
  rmst_oracle_basename     = rmst_oracle_basename,
  rmst_oracle_alt_basename = rmst_oracle_alt_basename,
  time_for_match_predben   = "P5",
  output_dir               = output_dir,
  rate_R                   = RATE_BOOTSTRAP_R
)

message(sprintf("Done — results in %s/resAll_CI_match%s.csv", output_dir, MATCH_TYPE))

