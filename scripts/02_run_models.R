# =============================================================================
# scripts/02_run_models.R
#
# Fit one R-based model on simulated data from setting 1 or 2.
#
# Available models:
#   "ALASSO"            — adaptive LASSO Cox proportional hazards (benchmark)
#   "InteractionForest" — Interaction Forest
#   "CSF"               — Causal Survival Forest
#
# For CoxCC / CoxTime use scripts/run_nn.py instead.
#
# Output per replication j — unified across all models and both settings:
#   pred_horizons.csv      — factual survival probabilities at P25, P5, P75
#   pred_horizons_alt.csv  — counterfactual survival probabilities
#   rmst_horizons.csv      — factual RMST at P25, P5, P75
#   rmst_horizons_alt.csv  — counterfactual RMST
#   benefit_surv.csv       — CSF only: patient-level CATE (survival scale)
#   benefit_rmst.csv       — CSF only: patient-level CATE (RMST scale)
#   pred_horizons_sf.csv   — CSF only: S-learner factual survival
#   rmst_horizons_sf.csv   — CSF only: S-learner factual RMST
#   rate_autoc.csv         — CSF only: AUTOC per horizon and target
#   pred.csv               — ALASSO only: full survival grid
#   selecVar.csv           — ALASSO only: selected variables
# =============================================================================

library(here)
source(here("config.R"))

# =============================================================================
# USER SETTINGS — edit here
# =============================================================================

SETTING <- 2          # 1 = Weibull/Biospear  |  2 = AFT/Friedman
MODEL   <- "CSF"   # "ALASSO" | "InteractionForest" | "CSF"
N_GRID  <- 200        # fine-grid resolution for ALASSO

# =============================================================================
# Paths
# =============================================================================
if (SETTING == 1) {
  data_root <- list.dirs(here("data", "sim_setting1"), recursive = FALSE)[1]
  res_root  <- here("results", "sim_setting1", MODEL)
} else if (SETTING == 2) {
  data_root <- list.dirs(here("data", "sim_setting2"), recursive = FALSE)[1]
  res_root  <- here("results", "sim_setting2", MODEL)
} else {
  stop("SETTING must be 1 or 2.")
}

sim_details <- basename(data_root)
dir.create(file.path(res_root, sim_details), recursive = TRUE, showWarnings = FALSE)
message(sprintf("Setting %d | model: %s | scenario: %s", SETTING, MODEL, sim_details))

# =============================================================================
# Load model helpers
# =============================================================================
if (MODEL == "ALASSO") {
  source(here("src", "models", "alasso.R"))
} else if (MODEL == "InteractionForest") {
  source(here("src", "models", "InteractionForest.R"))
} else if (MODEL == "CSF") {
  source(here("src", "models", "csf.R"))
} else {
  stop("MODEL must be 'ALASSO', 'InteractionForest', or 'CSF'.")
}

# =============================================================================
# Helper — read the per-replication horizon file
# =============================================================================
read_times <- function(j) {
  times <- read.csv(file.path(data_root, "surv_th",
                              paste0("simdata_", j), "times.csv"))
  list(
    P25      = as.numeric(times[["P25"]]),
    P5       = as.numeric(times[["P5"]]),
    P75      = as.numeric(times[["P75"]]),
    horizons = c(as.numeric(times[["P25"]]),
                 as.numeric(times[["P5"]]),
                 as.numeric(times[["P75"]]))
  )
}

# =============================================================================
# Main loop
# =============================================================================
for (j in ITERS) {
  
  message(sprintf("  replication %d / %d", j, max(ITERS)))
  
  dir_res <- file.path(res_root, sim_details, paste0("simdata_", j))
  dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)
  
  df_train <- read.csv(file.path(data_root, "all", paste0("simdata_", j), "train.csv"))
  df_test  <- read.csv(file.path(data_root, "all", paste0("simdata_", j), "test.csv"))
  t        <- read_times(j)
  
  # ===========================================================================
  # ALASSO
  # ===========================================================================
  if (MODEL == "ALASSO") {
    
    fine_grid <- seq(from = t$P75 / N_GRID, to = t$P75, length.out = N_GRID)
    
    model <- fit_alasso_coxph(df_train)
    write.csv(summary(model$resBM), file.path(dir_res, "selecVar.csv"))
    
    # factual
    newd          <- prepare_newdata(df_test, model$inter_cols)
    surv_full     <- predict_surv_grid(model$fit, newd, fine_grid)
    colnames(surv_full) <- as.character(fine_grid)
    surv_horizons <- predict_surv_grid(model$fit, newd, t$horizons)
    colnames(surv_horizons) <- TIME_LABELS
    rmst_horizons <- rmst_grid(surv_full, fine_grid, t$horizons)
    
    write.csv(as.data.frame(surv_full),     file.path(dir_res, "pred.csv"),          row.names = FALSE)
    write.csv(as.data.frame(surv_horizons), file.path(dir_res, "pred_horizons.csv"), row.names = FALSE)
    write.csv(as.data.frame(rmst_horizons), file.path(dir_res, "rmst_horizons.csv"), row.names = FALSE)
    
    # counterfactual
    df_test_alt     <- df_test; df_test_alt$Trt <- -df_test_alt$Trt
    newd_alt              <- prepare_newdata(df_test_alt, model$inter_cols)
    surv_alt_full         <- predict_surv_grid(model$fit, newd_alt, fine_grid)
    colnames(surv_alt_full) <- as.character(fine_grid)
    surv_alt_horizons     <- predict_surv_grid(model$fit, newd_alt, t$horizons)
    colnames(surv_alt_horizons) <- TIME_LABELS
    rmst_alt_horizons     <- rmst_grid(surv_alt_full, fine_grid, t$horizons)
    
    write.csv(as.data.frame(surv_alt_full),     file.path(dir_res, "pred_alt.csv"),          row.names = FALSE)
    write.csv(as.data.frame(surv_alt_horizons), file.path(dir_res, "pred_horizons_alt.csv"), row.names = FALSE)
    write.csv(as.data.frame(rmst_alt_horizons), file.path(dir_res, "rmst_horizons_alt.csv"), row.names = FALSE)
    
    gc()
    
    # ===========================================================================
    # Interaction Forest
    # ===========================================================================
  } else if (MODEL == "InteractionForest") {
    
    model <- fit_interaction_forest(df_train)
    rm(df_train); gc()
    
    df_test     <- subset(df_test, select = -X)
    df_test_alt <- df_test; df_test_alt$Trt <- -df_test_alt$Trt
    
    # factual
    res           <- predict_interaction_forest(model, df_test)
    surv          <- res$surv
    grid_times    <- res$times
    surv_horizons <- extract_at_horizons(as.matrix(surv), grid_times, t$horizons)
    colnames(surv_horizons) <- TIME_LABELS
    rmst_horizons <- rmst_grid(as.matrix(surv), grid_times, t$horizons)
    
    write.csv(surv,                           file.path(dir_res, "pred.csv"),          row.names = FALSE)
    write.csv(as.data.frame(surv_horizons),   file.path(dir_res, "pred_horizons.csv"), row.names = FALSE)
    write.csv(as.data.frame(rmst_horizons),   file.path(dir_res, "rmst_horizons.csv"), row.names = FALSE)
    
    # counterfactual
    res_alt           <- predict_interaction_forest(model, df_test_alt)
    surv_alt          <- res_alt$surv
    surv_alt_horizons <- extract_at_horizons(as.matrix(surv_alt), res_alt$times, t$horizons)
    colnames(surv_alt_horizons) <- TIME_LABELS
    rmst_alt_horizons <- rmst_grid(as.matrix(surv_alt), res_alt$times, t$horizons)
    
    write.csv(surv_alt,                          file.path(dir_res, "pred_alt.csv"),          row.names = FALSE)
    write.csv(as.data.frame(surv_alt_horizons),  file.path(dir_res, "pred_horizons_alt.csv"), row.names = FALSE)
    write.csv(as.data.frame(rmst_alt_horizons),  file.path(dir_res, "rmst_horizons_alt.csv"), row.names = FALSE)
    
    rm(model, df_test, df_test_alt, surv, surv_alt); gc()
    
    # ===========================================================================
    # Causal Survival Forest
    # ===========================================================================
  } else if (MODEL == "CSF") {
    
    drop_base   <- c("yy", "status", "X")
    X_train_csf <- as.matrix(df_train[, setdiff(names(df_train), c(drop_base, "Trt"))])
    X_test_csf  <- as.matrix(df_test[,  setdiff(names(df_test),  c(drop_base, "Trt"))])
    X_train_sf  <- as.matrix(df_train[, setdiff(names(df_train), drop_base)])
    X_test_sf   <- as.matrix(df_test[,  setdiff(names(df_test),  drop_base)])
    
    df_test_alt     <- df_test; df_test_alt$Trt <- -df_test_alt$Trt
    X_test_sf_alt   <- as.matrix(df_test_alt[, setdiff(names(df_test_alt), drop_base)])
    
    # Trt is {-0.5, 0.5}; shift to {0, 1} for grf
    Y_train <- df_train$yy;  W_train <- df_train$Trt + 0.5;  D_train <- as.integer(df_train$status)
    Y_test  <- df_test$yy;   W_test  <- df_test$Trt  + 0.5;  D_test  <- as.integer(df_test$status)
    
    sf_all <- fit_survival_forest(X_train_sf, Y_train, D_train)
    
    benefit_surv <- matrix(NA_real_, nrow = nrow(X_test_csf), ncol = 3, dimnames = list(NULL, TIME_LABELS))
    benefit_rmst <- benefit_surv
    surv_sf      <- matrix(NA_real_, nrow = nrow(X_test_sf),  ncol = 3, dimnames = list(NULL, TIME_LABELS))
    surv_sf_alt  <- surv_sf; rmst_sf <- surv_sf; rmst_sf_alt <- surv_sf
    rate_rows    <- list()
    
    for (k in seq_along(t$horizons)) {
      t_k  <- t$horizons[[k]]
      h_nm <- TIME_LABELS[[k]]
      
      surv_sf[, k]     <- surv_at_horizon(sf_all, X_test_sf,     t_k)
      surv_sf_alt[, k] <- surv_at_horizon(sf_all, X_test_sf_alt, t_k)
      rmst_sf[, k]     <- rmst_from_sf(sf_all,    X_test_sf,     t_k)
      rmst_sf_alt[, k] <- rmst_from_sf(sf_all,    X_test_sf_alt, t_k)
      
      csf_sp            <- fit_causal_survival_forest(X_train_csf, Y_train, W_train, D_train, t_k, "survival.probability")
      benefit_surv[, k] <- predict(csf_sp, X_test_csf)$predictions
      
      csf_rmst          <- fit_causal_survival_forest(X_train_csf, Y_train, W_train, D_train, t_k, "RMST")
      benefit_rmst[, k] <- predict(csf_rmst, X_test_csf)$predictions
      
      for (tgt in c("survival.probability", "RMST")) {
        score_forest <- if (tgt == "survival.probability") csf_sp else csf_rmst
        eval_forest  <- fit_causal_survival_forest(X_test_csf, Y_test, W_test, D_test, t_k, tgt)
        autoc_df     <- compute_autoc(score_forest, eval_forest, X_test_csf)
        rate_rows[[length(rate_rows) + 1L]] <- cbind(
          data.frame(target = tgt, horizon = h_nm, horizon_t = t_k, stringsAsFactors = FALSE),
          autoc_df
        )
      }
    }
    
    # CSF-specific outputs
    write.csv(as.data.frame(benefit_surv), file.path(dir_res, "benefit_surv.csv"),        row.names = FALSE)
    write.csv(as.data.frame(benefit_rmst), file.path(dir_res, "benefit_rmst.csv"),        row.names = FALSE)
    write.csv(as.data.frame(surv_sf),      file.path(dir_res, "pred_horizons_sf.csv"),    row.names = FALSE)
    write.csv(as.data.frame(rmst_sf),      file.path(dir_res, "rmst_horizons_sf.csv"),    row.names = FALSE)
    write.csv(do.call(rbind, rate_rows),   file.path(dir_res, "rate_autoc.csv"),           row.names = FALSE)
    
    # Unified outputs (same names as other models)
    write.csv(as.data.frame(benefit_surv), file.path(dir_res, "pred_horizons.csv"),       row.names = FALSE)
    write.csv(as.data.frame(benefit_rmst), file.path(dir_res, "rmst_horizons.csv"),       row.names = FALSE)
    write.csv(as.data.frame(surv_sf_alt),  file.path(dir_res, "pred_horizons_alt.csv"),   row.names = FALSE)
    write.csv(as.data.frame(rmst_sf_alt),  file.path(dir_res, "rmst_horizons_alt.csv"),   row.names = FALSE)
    
  }
  
  message(sprintf("    done replication %d", j))
}

message(sprintf("Complete — results in results/sim_setting%d/%s/", SETTING, MODEL))

