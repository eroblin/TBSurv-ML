# =============================================================================
# src/metrics/discrimination.R
# Metrics that quantify how well the model rank-orders patients.
#
#   compute_uno_c        : Uno's C-statistic for outcomes (overall, treated, untreated)
#   compute_c_for_benefit: van Klaveren's C-for-benefit on matched pairs
#   compute_rate         : RATE (AUTOC and QINI) on the full test set
# =============================================================================

# expected packages: survival, survAUC, Hmisc, grf
# (caller is expected to library() these; we don't here)

# -----------------------------------------------------------------------------
# Uno's C on outcomes, overall and by treatment arm.
#
# Returns a named numeric vector with three values.
# -----------------------------------------------------------------------------
compute_uno_c <- function(yy_train, status_train,
                          yy_test,  status_test,
                          treat_test, pred_test, time) {
  
  Surv.train <- survival::Surv(yy_train, status_train)
  t.star     <- time_from_label(time)
  pred_vec   <- as.numeric(unlist(pred_test[time]))
  
  # overall ------------------------------------------------------------------
  Surv.test <- survival::Surv(yy_test, status_test)
  c_overall <- survAUC::UnoC(Surv.train, Surv.test, -pred_vec, t.star)
  
  # treated subset -----------------------------------------------------------
  is_trt    <- treat_test == 1 | treat_test == 1.5   # accept either coding
  Surv.trt  <- survival::Surv(yy_test[is_trt], status_test[is_trt])
  c_trt     <- survAUC::UnoC(Surv.train, Surv.trt, -pred_vec[is_trt], t.star)
  
  # untreated subset ---------------------------------------------------------
  Surv.untrt <- survival::Surv(yy_test[!is_trt], status_test[!is_trt])
  c_untrt    <- survAUC::UnoC(Surv.train, Surv.untrt, -pred_vec[!is_trt], t.star)
  
  c(c = c_overall, c_trt = c_trt, c_untrt = c_untrt)
}

# -----------------------------------------------------------------------------
# C-for-benefit on matched (treated, control) pairs.
#
# Inputs:
#   tau_predicted : pairwise predicted benefit (averaged within pair).
#   tau_observed  : pairwise observed benefit (+1, 0, -1) using survival times
#                   and statuses of the two members of the pair.
#
# We compute tau_observed here using the original definition from
# van Klaveren (2018) as implemented in your legacy code.
# -----------------------------------------------------------------------------
compute_c_for_benefit <- function(tau_predicted, tau_observed) {
  Hmisc::rcorr.cens(tau_predicted, tau_observed)[[1]]
}

#' Helper: build pairwise observed benefit from a matched-pair dataset.
#'
#' Returns +1 if the control member of the pair had an event before the
#' treated member, -1 if the reverse, 0 if neither comparison can be made
#' (both censored, or censoring obscures the order).
build_tau_observed <- function(data, treated_idx, control_idx,
                               time_var = "yy", status_var = "status") {
  ctrl_event_first <- data[control_idx, status_var] == 1 &
    data[control_idx, time_var]   < data[treated_idx, time_var]
  trt_event_first  <- data[treated_idx, status_var] == 1 &
    data[treated_idx, time_var]   < data[control_idx, time_var]
  as.numeric(ctrl_event_first) - as.numeric(trt_event_first)
}


build_tau_observed_rmst <- function(data, treated_idx, control_idx,
                                    time_var = "yy", t_star) {
  yy_T <- pmin(data[treated_idx, time_var], t_star)
  yy_C <- pmin(data[control_idx, time_var], t_star)
  yy_C - yy_T
}
# -----------------------------------------------------------------------------
# RATE: Rank-Average Treatment Effect (Yadlowsky et al. 2023).
#
# Computes both AUTOC and QINI variants. Uses doubly-robust IPW scores valid
# under randomization. Priorities are the patient-level predicted benefit
# (already corrected for the within-pair sign flip).
#
# Input:
#   trt01      : 0/1 treatment indicator (NOT shifted)
#   status     : event indicator (used as outcome here, matching your code)
#   priorities : patient-level signed predicted benefit
# -----------------------------------------------------------------------------
compute_rate <- function(trt01, status, priorities, R = 200) {
  
  W     <- as.numeric(trt01 == 1)
  W.hat <- mean(W)
  Y     <- status
  
  DR.scores <- ifelse(W == 1, Y / W.hat, -Y / (1 - W.hat))
  
  rate_autoc <- tryCatch(
    grf::rank_average_treatment_effect.fit(
      DR.scores  = DR.scores,
      priorities = priorities,
      target     = "AUTOC",
      R          = R
    ),
    error = function(e) { cat("AUTOC error:", conditionMessage(e), "\n"); NULL }
  )
  
  rate_qini <- tryCatch(
    grf::rank_average_treatment_effect.fit(
      DR.scores  = DR.scores,
      priorities = priorities,
      target     = "QINI",
      R          = R
    ),
    error = function(e) { cat("QINI error:", conditionMessage(e), "\n"); NULL }
  )
  
  autoc_est <- if (!is.null(rate_autoc)) rate_autoc$estimate else NA_real_
  autoc_se  <- if (!is.null(rate_autoc)) rate_autoc$std.err  else NA_real_
  qini_est  <- if (!is.null(rate_qini))  rate_qini$estimate  else NA_real_
  qini_se   <- if (!is.null(rate_qini))  rate_qini$std.err   else NA_real_
  
  
  c(RATE_AUTOC = autoc_est, RATE_QINI = qini_est)
}