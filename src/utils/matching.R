# =============================================================================
# src/utils/matching.R
# Two matching strategies used in the simulations:
#   - "var"    : Mahalanobis distance on baseline covariates (MatchIt)
#   - "predben": rank-based matching on predicted benefit
# Both return a list with `treated` and `control` index vectors of equal length.
# =============================================================================
#' Mahalanobis matching on covariates.
#'
#' @param df_test    
#' @param treat_var  
#' @param drop_cols  
match_mahalanobis <- function(df_test,
                              treat_var = "Trt",
                              drop_cols = c("X", "yy", "status",
                                            "Trt", "pred_ben", "pred_benth")) {
  list_var <- setdiff(colnames(df_test), drop_cols)
  fml      <- as.formula(paste(treat_var, "~", paste(list_var, collapse = " + ")))
  
  fit <- MatchIt::matchit(
    formula  = fml,
    data     = df_test,
    distance = "mahalanobis",
    replace  = FALSE,
    estimand = NULL
  )
  m <- fit$match.matrix
  list(
    treated = as.numeric(rownames(m)),
    control = as.numeric(m[, 1])
  )
}

#' Rank-match on predicted benefit.
#'
#' Sort treated patients by their predicted benefit, sort controls by theirs,
#' then pair the i-th of each. Truncate to the smaller arm so pairs are 1:1.
#'
#' @param df_test  test set; must contain a column `pred_ben` already
#'                 (patient-level signed predicted benefit).
#' @param treat_var name of the treatment indicator column.
#' @param treated_value value indicating "treated" in `treat_var`
#'                      (1 if you kept original coding, 1.5 if shifted).
match_predben <- function(df_test,
                          treat_var     = "Trt",
                          treated_value = 1) {
  is_trt <- df_test[[treat_var]] == treated_value
  ind_A  <- which(is_trt)        # treated indices
  ind_B  <- which(!is_trt)       # control indices
  
  n_pairs <- min(length(ind_A), length(ind_B))
  
  # order each arm by predicted benefit (ascending)
  ord_A <- order(df_test[ind_A, "pred_ben"])
  ord_B <- order(df_test[ind_B, "pred_ben"])
  
  list(
    treated = ind_A[ord_A][seq_len(n_pairs)],
    control = ind_B[ord_B][seq_len(n_pairs)]
  )
}

#' Dispatcher.
make_matches <- function(df_test, match_type, treat_var = "Trt",
                         treated_value = 1, drop_cols = NULL) {
  if (match_type == "var") {
    if (is.null(drop_cols)) {
      drop_cols <- c("X", "yy", "status", "Trt", "pred_ben", "pred_benth")
    }
    return(match_mahalanobis(df_test, treat_var, drop_cols))
  }
  if (match_type == "predben") {
    return(match_predben(df_test, treat_var, treated_value))
  }
  stop("Unknown match_type: ", match_type)
}