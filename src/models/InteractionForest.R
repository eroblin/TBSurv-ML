# =============================================================================
# src/models/interaction_forest.R
# Interaction forest model for survival data.
# =============================================================================

library(diversityForest)

# -----------------------------------------------------------------------------
# Extract survival probabilities at given horizons from a full survival matrix
# -----------------------------------------------------------------------------
extract_at_horizons <- function(surv_full, grid_times, horizons) {
  out <- matrix(1, nrow = nrow(surv_full), ncol = length(horizons))
  for (k in seq_along(horizons)) {
    idx <- which(grid_times <= horizons[k])
    if (length(idx) > 0) out[, k] <- surv_full[, max(idx)]
  }
  out
}

# -----------------------------------------------------------------------------
# RMST at given horizons via trapezoidal rule
# -----------------------------------------------------------------------------
rmst_grid <- function(surv_full, grid_times, horizons) {
  times  <- c(0, as.numeric(grid_times))
  S      <- cbind(1, surv_full)
  n_subj <- nrow(S)
  
  out <- matrix(0, nrow = n_subj, ncol = length(horizons))
  for (k in seq_along(horizons)) {
    tstar <- horizons[k]
    mask  <- times < tstar
    t_sub <- c(times[mask], tstar)
    
    if (tstar <= min(times)) {
      s_at_tstar <- rep(1, n_subj)
    } else if (tstar >= max(times)) {
      s_at_tstar <- S[, ncol(S)]
    } else {
      i_hi       <- which(times >= tstar)[1]
      i_lo       <- i_hi - 1
      w          <- (tstar - times[i_lo]) / (times[i_hi] - times[i_lo])
      s_at_tstar <- (1 - w) * S[, i_lo] + w * S[, i_hi]
    }
    
    S_sub    <- cbind(S[, mask, drop = FALSE], s_at_tstar)
    dt       <- diff(t_sub)
    avgS     <- (S_sub[, -1, drop = FALSE] + S_sub[, -ncol(S_sub), drop = FALSE]) / 2
    out[, k] <- as.numeric(avgS %*% dt)
  }
  colnames(out) <- c("P25", "P5", "P75")
  out
}

# -----------------------------------------------------------------------------
# Fit interaction forest
# -----------------------------------------------------------------------------
fit_interaction_forest <- function(df_train, max_depth = 5) {
  df_train$Trt <- as.factor(df_train$Trt)
  df_train     <- subset(df_train, select = -c(X))
  
  interactionfor(dependent.variable.name = "yy",
                 data                    = df_train,
                 status.variable.name    = "status",
                 importance              = "none",
                 simplify.large.n        = TRUE,
                 max.depth               = max_depth,
                 verbose                 = TRUE)
}

# -----------------------------------------------------------------------------
# Predict survival on test data
# -----------------------------------------------------------------------------
predict_interaction_forest <- function(model, df_test) {
  df_test$Trt    <- as.factor(df_test$Trt)
  pred           <- predict(model, data = df_test, num.threads = 24)
  surv           <- as.data.frame(pred$survival)
  colnames(surv) <- pred$unique.death.times
  list(surv = surv, times = pred$unique.death.times)
}