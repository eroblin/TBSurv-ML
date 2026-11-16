# =============================================================================
# output_results.R
# Per-iteration orchestrator.
#
# Computes:
#   - survival-probability discrimination + calibration
#   - RMST-based discrimination + calibration ('_rmst' suffix)
#   - PB family (both scales) on the full test set
#
# CSF path (auto-detected via the global MODEL == "CSF"):
#   - patient-level signed benefit is read directly from
#       <dir_plot>/benefit_surv.csv   (survival scale)
#       <dir_plot>/benefit_rmst.csv   (RMST scale)
#   - pred_alt / rmst_alt are NOT used (they're only needed to construct the
#     benefit via (p - p_alt) * (-1)^(1-Trt) for the non-CSF path).
# =============================================================================

output_results <- function(treated, control,
                           df_train, df_test,
                           pred, pred_alt,
                           pth,  pth_alt,
                           rmst, rmst_alt,
                           rmst_th, rmst_th_alt,
                           time,
                           dir_plot   = NULL,
                           treat_var  = "Trt",
                           rate_R     = 200,
                           cell_ok    = TRUE) {
  
  trt01 <- df_test[[treat_var]]
  
  is_csf <- exists("MODEL", envir = .GlobalEnv) &&
    identical(get("MODEL", envir = .GlobalEnv), "CSF")
  
  # ---- 0. patient-level signed predicted benefit (both scales) -----------
  if (is_csf) {
    benefit_surv <- read.csv(file.path(dir_plot, "benefit_surv.csv"))
    benefit_rmst <- read.csv(file.path(dir_plot, "benefit_rmst.csv"))
    pred_ben_vec <- as.numeric(benefit_surv[[time]])
    rmst_ben_vec <- as.numeric(benefit_rmst[[time]])
  } else {
    pred_ben_vec <- patient_pred_ben(pred, pred_alt, time, trt01)
    rmst_ben_vec <- patient_pred_ben(rmst, rmst_alt, time, trt01)
  }

  # ────────────────────────────────────────────────────────────────────────
  
  # Oracle benefit is always built from the oracle pair.
  predth_vec      <- patient_pred_ben(pth,     pth_alt,     time, trt01)
  rmst_ben_th_vec <- patient_pred_ben(rmst_th, rmst_th_alt, time, trt01)
  
  # ────────────────────────────────────────────────────────────────────────
  
  # ---- 1. pairwise predicted/observed benefit on matched pairs ----------
  tau_predicted <- pairwise_pred_ben(pred_ben_vec, treated, control)
  
  cleaned <- drop_na_pairs(tau_predicted, treated, control)
  tau_predicted <- cleaned$tau
  treated       <- cleaned$treated
  control       <- cleaned$control
  

  
  tau_predicted_oracle <- pairwise_pred_ben(predth_vec, treated, control)
  tau_observed         <- build_tau_observed(df_test, treated, control)
  
  tau_predicted_rmst        <- pairwise_pred_ben(rmst_ben_vec,    treated, control)
  tau_predicted_rmst_oracle <- pairwise_pred_ben(rmst_ben_th_vec, treated, control)
  tau_observed_rmst         <- build_tau_observed_rmst(df_test, treated, control,
                                                       t_star = time_from_label(time))

  # ---- 2. discrimination ------------------------------------------------
  uno_c <- compute_uno_c(
    yy_train     = df_train$yy,
    status_train = df_train$status,
    yy_test      = df_test$yy,
    status_test  = df_test$status,
    treat_test   = df_test[[treat_var]],
    pred_test    = pred,
    time         = time
  )
  
  c_ben         <- compute_c_for_benefit(tau_predicted, tau_observed)
  c_ben_th      <- tryCatch(compute_c_for_benefit(tau_predicted_oracle, tau_observed),
                            error = function(e) NA_real_)
  c_ben_rmst    <- tryCatch(compute_c_for_benefit(tau_predicted_rmst, tau_observed),
                            error = function(e) NA_real_)
  c_ben_rmst_th <- tryCatch(compute_c_for_benefit(tau_predicted_rmst_oracle, tau_observed),
                            error = function(e) NA_real_)
  
  rate_vals <- compute_rate(
    trt01      = trt01,
    status     = df_test$status,
    priorities = pred_ben_vec,
    R          = rate_R
  )
  
  rate_vals_rmst <- compute_rate(
    trt01      = trt01,
    status     = df_test$status,
    priorities = rmst_ben_vec,
    R          = rate_R
  )
  
  # ---- 3. calibration ---------------------------------------------------
  na_cal <- list(Eavg = NA_real_, E50 = NA_real_, tau_smoothed = NA_real_)
  
  cal <- tryCatch(
    compute_calibration_for_benefit(tau_predicted, tau_observed),
    error = function(e) na_cal)
  cal_th <- tryCatch(
    compute_calibration_for_benefit(tau_predicted_oracle, tau_observed),
    error = function(e) na_cal)
  rmse <- tryCatch(
    compute_rmse_vs_oracle(tau_predicted, tau_predicted_oracle),
    error = function(e) NA_real_)
  
  cal_rmst <- tryCatch(
    compute_calibration_for_benefit(tau_predicted_rmst, tau_observed_rmst),
    error = function(e) na_cal)
  cal_rmst_th <- tryCatch(
    compute_calibration_for_benefit(tau_predicted_rmst_oracle, tau_observed_rmst),
    error = function(e) na_cal)
  rmse_rmst <- tryCatch(
    compute_rmse_vs_oracle(tau_predicted_rmst, tau_predicted_rmst_oracle),
    error = function(e) NA_real_)
  
  if (!is.null(dir_plot) && time == "P5" && !all(is.na(cal$tau_smoothed))) {
    p <- plot_calibration_for_benefit(tau_predicted, cal$tau_smoothed)
    save_calibration_plot(p, dir_plot)
  }
  
  # ---- 4. PB family on full test set ------------------------------------
  if (cell_ok) {
    pb_vals <- pb_all(
      pred_ben = pred_ben_vec,
      trt01    = trt01,
      time     = df_test$yy,
      status   = df_test$status,
      t.star   = time_from_label(time)
    )
  } else {
    pb_names <- c(
      "PB.est","PB.se","PB.lo","PB.hi",
      "PB_MvT.est","PB_MvT.se","PB_MvT.lo","PB_MvT.hi",
      "PB_MvC.est","PB_MvC.se","PB_MvC.lo","PB_MvC.hi",
      "PB_treat_all.est","PB_treat_all.se","PB_treat_all.lo","PB_treat_all.hi",
      "PB_RMST.est","PB_RMST.se","PB_RMST.lo","PB_RMST.hi",
      "PB_RMST_MvT.est","PB_RMST_MvT.se","PB_RMST_MvT.lo","PB_RMST_MvT.hi",
      "PB_RMST_MvC.est","PB_RMST_MvC.se","PB_RMST_MvC.lo","PB_RMST_MvC.hi",
      "PB_RMST_treat_all.est","PB_RMST_treat_all.se",
      "PB_RMST_treat_all.lo","PB_RMST_treat_all.hi"
    )
    pb_vals <- setNames(rep(NA_real_, length(pb_names)), pb_names)
  }
  
  
  # ---- 5. assemble ------------------------------------------------------
  legacy <- c(
    c                  = unname(uno_c["c"]),
    c_trt              = unname(uno_c["c_trt"]),
    c_untrt            = unname(uno_c["c_untrt"]),
    c_ben              = c_ben,
    Eavg.for.benefit   = cal$Eavg,
    E50.for.benefit    = cal$E50,
    c_benth            = c_ben_th,
    Eavgth.for.benefit = cal_th$Eavg,
    E50th.for.benefit  = cal_th$E50,
    RMSE               = rmse,
    RATE_AUTOC         = unname(rate_vals["RATE_AUTOC"]),
    RATE_QINI          = unname(rate_vals["RATE_QINI"])
  )
  
  rmst_block <- c(
    c_ben_rmst              = c_ben_rmst,
    Eavg.for.benefit_rmst   = cal_rmst$Eavg,
    E50.for.benefit_rmst    = cal_rmst$E50,
    c_benth_rmst            = c_ben_rmst_th,
    Eavgth.for.benefit_rmst = cal_rmst_th$Eavg,
    E50th.for.benefit_rmst  = cal_rmst_th$E50,
    RMSE_rmst               = rmse_rmst,
    RATE_AUTOC_rmst         = unname(rate_vals_rmst["RATE_AUTOC"]),
    RATE_QINI_rmst          = unname(rate_vals_rmst["RATE_QINI"])
  )
  
  c(legacy, rmst_block, pb_vals)
}