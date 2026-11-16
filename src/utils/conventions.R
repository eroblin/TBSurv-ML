# =============================================================================
# src/utils/conventions.R
# =============================================================================

# ---- Conventions ------------------------------------------------------------
#
# 1. Treatment coding inside df_test:
#       Raw CSVs have Trt in {-0.5, +0.5}.
#       The pipeline shifts to {0, 1} via `df_test$Trt <- df_test$Trt + 0.5`.
#       After the shift, every metric function uses Trt directly as 0/1.
#
# 2. Predicted benefit, two different uses:
#
#    a) PATIENT-LEVEL predicted benefit (used by PB and RATE):
#         pred_ben_i = pred[time]_i  -  pred_alt[time]_i
#       where pred is "predicted survival under the regimen the patient
#       was randomised to" and pred_alt is "under the counterfactual regimen".
#       Sign convention here: positive => model thinks treatment helps patient i.
#       This needs a sign flip for control patients because for them
#       pred = S(control) and pred_alt = S(treatment), so the raw difference
#       has the wrong sign. The flip is `(-1)^(1 - trt01)`:
#         trt01 == 1 -> multiplier  +1
#         trt01 == 0 -> multiplier  -1
#
#    b) PAIRWISE predicted benefit (used by C-for-benefit, calibration loess):
#       For a matched (treated, control) pair the pairwise predicted benefit is
#         tau_pred = (pred_ben_treated + pred_ben_control) / 2
#       This is what `output_results` computed before; we keep it in the
#       discrimination/calibration files.
#
# 3. Time-horizon column names "P25", "P5", "P75"
# -----------------------------------------------------------------------------

#' Parse the numeric time horizon from a column label like "P25" -> 25.
time_from_label <- function(label) {
  as.numeric(sub("^[A-Za-z]+", "", label))
}
#' Patient-level predicted benefit (positive = model thinks treatment helps).
#'
#' @param pred       data.frame of predicted survival probabilities at each
#'                   horizon under the *factual* regimen.
#' @param pred_alt   same shape, but under the *counterfactual* regimen.
#' @param time       column label, e.g. "P5".
#' @param trt01      0/1 treatment indicator (NOT the shifted version).
#' @return numeric vector of length nrow(pred).
patient_pred_ben <- function(pred, pred_alt, time, trt01) {
  raw <- as.numeric(unlist(pred[time] - pred_alt[time]))
  # flip sign for controls so positive always means "treatment helps"
  raw * ((-1) ^ (1 - trt01))
}

#' Pairwise predicted benefit averaged over a (treated, control) pair.
pairwise_pred_ben <- function(pred_ben_vec, treated_idx, control_idx) {
  (pred_ben_vec[treated_idx] + pred_ben_vec[control_idx]) / 2
}

#' Drop NA pairs from a parallel set of vectors.
#' Returns a list with the cleaned vectors and the count dropped.
drop_na_pairs <- function(tau, treated, control) {
  bad <- which(is.na(tau))
  if (length(bad) == 0) {
    return(list(tau = tau, treated = treated, control = control, n_dropped = 0))
  }
  list(
    tau       = tau[-bad],
    treated   = treated[-bad],
    control   = control[-bad],
    n_dropped = length(bad)
  )
}

#' Wrap a point estimate and SE into est / se / 95% CI.
with_ci <- function(est, se) {
  c(est = unname(est),
    se  = unname(se),
    lo  = unname(est - 1.96 * se),
    hi  = unname(est + 1.96 * se))
}