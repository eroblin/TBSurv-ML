# =============================================================================
# src/reporting/summarise.R
# Aggregate per-iteration results
# =============================================================================

#' Mean and SD per metric column across simulations.
#'
#' @param mat          numeric matrix, rows = simulation iterations,
#'                     cols = metrics.
#' @param metric_names column names to attach.
compute_summary <- function(mat, metric_names) {
  mn <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, sd, na.rm = TRUE)
  
  out <- rbind(mn, sdv)
  rownames(out) <- c("Mean", "SD")
  colnames(out) <- metric_names
  out
}

#' Combine summaries for the three time horizons into one data.frame.
format_output <- function(res25, res5, res75, metric_names) {
  df25 <- as.data.frame(t(res25))
  colnames(df25) <- paste0("P25_", c("Mean", "SD"))
  
  df5 <- as.data.frame(t(res5))
  colnames(df5) <- paste0("P5_", c("Mean", "SD"))
  
  df75 <- as.data.frame(t(res75))
  colnames(df75) <- paste0("P75_", c("Mean", "SD"))
  
  out <- cbind(df25, df5, df75)
  round(out, 3)
}