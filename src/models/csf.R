# =============================================================================
# src/models/csf.R
# Causal survival forest and survival forest helpers.
# =============================================================================

library(grf)

# -----------------------------------------------------------------------------
# Survival probability S(horizon | X) from a survival_forest, per subject
# -----------------------------------------------------------------------------
surv_at_horizon <- function(sf_fit, X_new, horizon) {
  pr     <- predict(sf_fit, X_new)
  S      <- pr$predictions
  t_grid <- pr$failure.times
  idx    <- sum(t_grid <= horizon)
  if (idx == 0L) rep(1, nrow(X_new)) else S[, idx]
}

# -----------------------------------------------------------------------------
# RMST at horizon from a survival_forest via trapezoidal rule
# -----------------------------------------------------------------------------
rmst_from_sf <- function(sf_fit, X_new, horizon) {
  pr     <- predict(sf_fit, X_new)
  S      <- pr$predictions
  t_grid <- pr$failure.times
  
  keep         <- t_grid <= horizon
  S_at_0       <- rep(1, nrow(S))
  S_kept       <- S[, keep, drop = FALSE]
  S_at_horizon <- if (any(keep)) S_kept[, ncol(S_kept)] else S_at_0
  
  t_use <- c(0, t_grid[keep], horizon)
  S_use <- cbind(S_at_0, S_kept, S_at_horizon)
  dt    <- diff(t_use)
  S_avg <- (S_use[, -1, drop = FALSE] + S_use[, -ncol(S_use), drop = FALSE]) / 2
  as.numeric(S_avg %*% dt)
}

# -----------------------------------------------------------------------------
# Fit survival forest (S-learner, Trt included in X)
# -----------------------------------------------------------------------------
fit_survival_forest <- function(X_train, Y_train, D_train,
                                num_trees = 2000, min_node_size = 5,
                                alpha = 0.05, num_threads = 24) {
  survival_forest(
    X             = X_train,
    Y             = Y_train,
    D             = D_train,
    num.trees     = num_trees,
    min.node.size = min_node_size,
    alpha         = alpha,
    honesty       = TRUE,
    num.threads   = num_threads
  )
}

# -----------------------------------------------------------------------------
# Fit causal survival forest (Trt excluded from X)
# -----------------------------------------------------------------------------
fit_causal_survival_forest <- function(X_train, Y_train, W_train, D_train,
                                       horizon, target = "survival.probability",
                                       num_trees = 2000, min_node_size = 5,
                                       alpha = 0.05, num_threads = 24) {
  causal_survival_forest(
    X               = X_train,
    Y               = Y_train,
    W               = W_train,
    D               = D_train,
    W.hat           = mean(W_train),
    target          = target,
    horizon         = horizon,
    num.trees       = num_trees,
    min.node.size   = min_node_size,
    alpha           = alpha,
    honesty         = TRUE,
    compute.oob.predictions = TRUE,
    num.threads     = num_threads
  )
}

# -----------------------------------------------------------------------------
# Compute RATE/AUTOC given a scoring forest and an evaluation forest
# -----------------------------------------------------------------------------
compute_autoc <- function(score_forest, eval_forest, X_test) {
  priorities <- predict(score_forest, X_test)$predictions
  rate_obj   <- rank_average_treatment_effect(eval_forest, priorities, target = "AUTOC")
  data.frame(
    autoc     = rate_obj$estimate,
    std_err   = rate_obj$std.err,
    ci_lower  = rate_obj$estimate - 1.96 * rate_obj$std.err,
    ci_upper  = rate_obj$estimate + 1.96 * rate_obj$std.err,
    stringsAsFactors = FALSE
  )
}