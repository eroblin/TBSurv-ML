# =============================================================================
# run_scenario.R
# Generic scenario runner. Iterates over simulation replications, loads
# train/test data + predictions + oracle (survival probability AND RMST),
# runs matching, calls output_results() for three time horizons, and writes
# summary CSVs.
#
# When the four (Trt, pred_ben sign) cells are not all populated above
# `min_cell_predben`, PB-related metric columns are set to NA for that
# iteration but C / Eavg / RMSE / RATE columns are computed normally.
#
# CSF: when the global MODEL == "CSF", pred_alt / rmst_alt are not loaded
# (they don't exist on disk for CSF) and the patient-level signed benefit
# is read from benefit_surv.csv / benefit_rmst.csv. output_results() does
# the same detection.
# =============================================================================
main_dir <- "/Users/elr4013/Desktop/ISMB/"
setwd(main_dir)
source(here::here("src", "utils", "matching.R"))
source(here::here("src", "utils", "conventions.R"))
source(here::here("src", "metrics", "calibration.R"))
source(here::here("src", "metrics", "discrimination.R"))
source(here::here("src", "metrics", "decision_accuracy.R"))
source(here::here("src", "metrics", "output_results.R"))
source(here::here("src", "reporting", "summarise.R"))

run_scenario <- function(scenario_name,
                         data_root,
                         pred_root,
                         oracle_root,
                         iters,
                         match_type,
                         pred_basename          = "pred",
                         pred_alt_basename      = "pred_alt",
                         oracle_basename        = "predth_times",
                         oracle_alt_basename    = "predth_alt_times",
                         rmst_basename          = "rmst_horizons",
                         rmst_alt_basename      = "rmst_horizons_alt",
                         rmst_oracle_basename   = "rmst_test_times",
                         rmst_oracle_alt_basename = "rmst_test_alt_times",
                         time_for_match_predben = "P5",
                         output_dir             = pred_root,
                         rate_R                 = 200,
                         min_cell_predben       = 5) {
  
  is_csf <- exists("MODEL", envir = .GlobalEnv) &&
    identical(get("MODEL", envir = .GlobalEnv), "CSF")
  if (is_csf) {
    cat("[run_scenario] CSF detected (MODEL=='CSF'): pred_alt/rmst_alt skipped, ",
        "benefit_*.csv used.\n", sep = "")
  }
  
  res25_all <- NULL
  res5_all  <- NULL
  res75_all <- NULL
  skipped_iters_4groups <- integer(0)
  skipped_iters_error   <- integer(0)
  expected_names        <- NULL
  
  for (i in iters) {
    cat(sprintf("[%s, match=%s] iter %d\n", scenario_name, match_type, i))
    
    # ---- load data ------------------------------------------------------
    df_train <- read.csv(file.path(data_root, "all",
                                   paste0("simdata_", i), "train.csv"))
    df_test  <- read.csv(file.path(data_root, "all",
                                   paste0("simdata_", i), "test.csv"))
    df_test$Trt <- df_test$Trt + 0.5
    
    # survival-probability predictions
    p <- read.csv(file.path(pred_root, paste0("simdata_", i),
                            paste0(pred_basename, ".csv")))
    if (is_csf) {
      p_alt <- NULL
    } else {
      p_alt <- read.csv(file.path(pred_root, paste0("simdata_", i),
                                  paste0(pred_alt_basename, ".csv")))
    }
    
    # survival-probability oracle
    pth     <- read.csv(file.path(oracle_root, paste0("simdata_", i),
                                  paste0(oracle_basename, ".csv")))
    pth_alt <- read.csv(file.path(oracle_root, paste0("simdata_", i),
                                  paste0(oracle_alt_basename, ".csv")))
    
    # RMST predictions (model)
    rmst <- read.csv(file.path(pred_root, paste0("simdata_", i),
                               paste0(rmst_basename, ".csv")))
    if (is_csf) {
      rmst_alt <- NULL
    } else {
      rmst_alt <- read.csv(file.path(pred_root, paste0("simdata_", i),
                                     paste0(rmst_alt_basename, ".csv")))
    }
    
    # RMST oracle
    rmst_th     <- read.csv(file.path(oracle_root, paste0("simdata_", i),
                                      paste0(rmst_oracle_basename, ".csv")))
    rmst_th_alt <- read.csv(file.path(oracle_root, paste0("simdata_", i),
                                      paste0(rmst_oracle_alt_basename, ".csv")))
    
    # ---- patient-level pred_ben (survival scale) for the cell guard ----
    if (is_csf) {
      bs <- read.csv(file.path(pred_root, paste0("simdata_", i), "benefit_surv.csv"))
      df_test$pred_ben <- as.numeric(bs[[time_for_match_predben]])
    } else {
      df_test$pred_ben <- as.numeric(unlist(
        (p[time_for_match_predben] - p_alt[time_for_match_predben]) *
          ((-1) ^ (1 - df_test$Trt))
      ))
    }
    
    # ---- 4-group cell check (flag only; iter still runs) ---------------
    n11 <- sum(df_test$Trt == 1 & df_test$pred_ben > 0)
    n10 <- sum(df_test$Trt == 1 & df_test$pred_ben < 0)
    n01 <- sum(df_test$Trt == 0 & df_test$pred_ben > 0)
    n00 <- sum(df_test$Trt == 0 & df_test$pred_ben < 0)
    cell_ok <- min(n11, n10, n01, n00) >= min_cell_predben
    
    if (!cell_ok) {
      cat(sprintf("  >> iter %d: cells (T1+,T1-,T0+,T0-)=(%d,%d,%d,%d) < %d, PB cols -> NA\n",
                  i, n11, n10, n01, n00, min_cell_predben))
      skipped_iters_4groups <- c(skipped_iters_4groups, i)
    }
    
    # ---- compute metrics with tryCatch ---------------------------------
    computed <- tryCatch({
      pairs <- make_matches(df_test, match_type,
                            treat_var     = "Trt",
                            treated_value = 1)
      treatedR <- pairs$treated
      controlR <- pairs$control
      
      dir_plot <- file.path(pred_root, paste0("simdata_", i))
      
      list(
        r25 = output_results(treatedR, controlR, df_train, df_test,
                             p, p_alt, pth, pth_alt,
                             rmst, rmst_alt, rmst_th, rmst_th_alt,
                             "P25", dir_plot, rate_R = rate_R,
                             cell_ok = cell_ok),
        r5  = output_results(treatedR, controlR, df_train, df_test,
                             p, p_alt, pth, pth_alt,
                             rmst, rmst_alt, rmst_th, rmst_th_alt,
                             "P5",  dir_plot, rate_R = rate_R,
                             cell_ok = cell_ok),
        r75 = output_results(treatedR, controlR, df_train, df_test,
                             p, p_alt, pth, pth_alt,
                             rmst, rmst_alt, rmst_th, rmst_th_alt,
                             "P75", dir_plot, rate_R = rate_R,
                             cell_ok = cell_ok)
      )
    },error = function(e) {
      cat(sprintf("  >> iter %d: error during matching/metrics, NA row\n", i))
      cat("     msg:", conditionMessage(e), "\n")
      NULL
    })
    
    res25 <- NULL; res5 <- NULL; res75 <- NULL
    
    if (is.null(computed)) {
      skipped_iters_error <- c(skipped_iters_error, i)
    } else {
      res25 <- computed$r25
      res5  <- computed$r5
      res75 <- computed$r75
      
      if (is.null(expected_names)) {
        expected_names <- names(res5)
      } else if (!identical(names(res5), expected_names)) {
        cat(sprintf("  >> name mismatch at iter %d\n", i))
        print(setdiff(names(res5), expected_names))
        print(setdiff(expected_names, names(res5)))
        stop("Aborting at first mismatch.")
      }
    }
    
    # ---- append column for this iter (real values or NA placeholder) ---
    if (is.null(res5)) {
      na_col <- if (!is.null(expected_names))
        rep(NA_real_, length(expected_names))
      else NA_real_
      res25_all <- cbind(res25_all, na_col)
      res5_all  <- cbind(res5_all,  na_col)
      res75_all <- cbind(res75_all, na_col)
    } else {
      res25_all <- cbind(res25_all, res25)
      res5_all  <- cbind(res5_all,  res5)
      res75_all <- cbind(res75_all, res75)
    }
  }
  
  if (is.null(expected_names)) {
    cat(sprintf("[%s, match=%s] every iter errored, nothing to write.\n",
                scenario_name, match_type))
    return(invisible(NULL))
  }
  
  # ---- backfill any short placeholder columns ------------------------
  k <- length(expected_names)
  fix_mat <- function(mat) {
    n_iter <- length(iters)
    if (is.null(mat) || ncol(mat) != n_iter) return(mat)
    out <- matrix(NA_real_, nrow = k, ncol = n_iter)
    for (j in seq_len(n_iter)) {
      col <- mat[, j]
      if (length(col) == k) out[, j] <- col
    }
    out
  }
  res25_all <- fix_mat(res25_all)
  res5_all  <- fix_mat(res5_all)
  res75_all <- fix_mat(res75_all)
  
  res25_all <- t(res25_all); colnames(res25_all) <- expected_names
  res5_all  <- t(res5_all);  colnames(res5_all)  <- expected_names
  res75_all <- t(res75_all); colnames(res75_all) <- expected_names
  rownames(res25_all) <- paste0("simdata_", iters)
  rownames(res5_all)  <- paste0("simdata_", iters)
  rownames(res75_all) <- paste0("simdata_", iters)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  write.csv(as.data.frame(res5_all),
            file.path(output_dir,
                      sprintf("res5All_match%s.csv", match_type)),
            row.names = TRUE)
  
  sum25 <- compute_summary(res25_all, expected_names)
  sum5  <- compute_summary(res5_all,  expected_names)
  sum75 <- compute_summary(res75_all, expected_names)
  
  resR_final <- format_output(sum25, sum5, sum75, expected_names)
  
  # ---- bottom rows: skip counts (broadcast across all columns) -------
  n_4g <- length(skipped_iters_4groups)
  n_er <- length(skipped_iters_error)
  skip_row_4g <- data.frame(matrix(n_4g, nrow = 1, ncol = ncol(resR_final)))
  skip_row_er <- data.frame(matrix(n_er, nrow = 1, ncol = ncol(resR_final)))
  colnames(skip_row_4g) <- colnames(resR_final)
  colnames(skip_row_er) <- colnames(resR_final)
  rownames(skip_row_4g) <- "n_PB_NA_no4groups"
  rownames(skip_row_er) <- "n_skipped_runtime_error"
  resR_final <- rbind(resR_final, skip_row_4g, skip_row_er)
  
  write.csv(resR_final,
            file.path(output_dir,
                      sprintf("resAll_CI_match%s.csv", match_type)),
            row.names = TRUE)
  
  n_used <- length(iters) - length(skipped_iters_error)
  cat(sprintf("[%s, match=%s] done. %d / %d iters computed; PB NA in %d; runtime errors in %d.\n",
              scenario_name, match_type,
              n_used, length(iters),
              n_4g, n_er))
  if (n_4g > 0) cat("  PB-NA iters    :", paste(skipped_iters_4groups, collapse = ", "), "\n")
  if (n_er > 0) cat("  error iters    :", paste(skipped_iters_error, collapse = ", "), "\n")
  
  invisible(resR_final)
}