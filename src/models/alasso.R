# =============================================================================
# src/models/alasso.R
# ALASSO-Cox model: fitting, prediction, and RMST utilities.
# =============================================================================

library(survival)
library(biospear)

# -----------------------------------------------------------------------------
# Predict survival on a time grid for a fitted Cox model
# -----------------------------------------------------------------------------
predict_surv_grid <- function(coxfit, newdata, times) {
  s   <- survfit(coxfit, newdata = newdata, se.fit = FALSE)
  out <- matrix(1, nrow = nrow(newdata), ncol = length(times))
  for (k in seq_along(times)) {
    idx <- which(s$time <= times[k])
    if (length(idx) > 0) out[, k] <- s$surv[max(idx), ]
  }
  out
}

# -----------------------------------------------------------------------------
# RMST at given horizons via trapezoidal rule
# -----------------------------------------------------------------------------
rmst_grid <- function(surv_grid, fine_grid, horizons) {
  times  <- c(0, as.numeric(fine_grid))
  S      <- cbind(1, surv_grid)
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
# Fit ALASSO-Cox via BMsel
# -----------------------------------------------------------------------------
fit_alasso_coxph <- function(df_train) {
  genes <- setdiff(colnames(df_train), c("X", "yy", "status", "Trt"))
  
  resBM <- BMsel(df_train,
                 x      = genes,
                 y      = c("yy", "status"),
                 tt     = "Trt",
                 inter  = TRUE,
                 std.x  = TRUE,
                 std.i  = TRUE,
                 std.tt = FALSE,
                 method = "alassoU",
                 folds  = 5,
                 trace  = FALSE)
  
  coef_tab <- summary(resBM, show = FALSE)
  if (is.list(coef_tab)) coef_tab <- coef_tab[["alassoU"]]
  if (is.null(coef_tab)) stop("Could not extract alassoU coefficients from BMsel output.")
  
  selected <- rownames(coef_tab)[coef_tab[, 1] != 0]
  selected <- setdiff(selected, "(Intercept)")
  if (length(selected) == 0) stop("No variables selected by ALASSO for this simulation.")
  
  d          <- df_train
  inter_cols <- grep(":", selected, value = TRUE)
  for (nm in inter_cols) {
    parts    <- strsplit(nm, ":", fixed = TRUE)[[1]]
    d[[nm]]  <- d[[parts[1]]] * d[[parts[2]]]
  }
  
  fct <- as.formula(paste0("Surv(yy, status) ~ ",
                           paste0("`", selected, "`", collapse = " + ")))
  fit <- coxph(fct, data = d, init = coef_tab[selected, 1], iter = 0)
  list(fit = fit, selected = selected, inter_cols = inter_cols, resBM = resBM)
}

# -----------------------------------------------------------------------------
# Prepare newdata by adding interaction columns
# -----------------------------------------------------------------------------
prepare_newdata <- function(df, inter_cols) {
  for (nm in inter_cols) {
    parts   <- strsplit(nm, ":", fixed = TRUE)[[1]]
    df[[nm]] <- df[[parts[1]]] * df[[parts[2]]]
  }
  df
}