# =============================================================================
# src/metrics/calibration.R
# Calibration-for-benefit on matched pairs, plus the calibration plot.
# =============================================================================
# Calibration-for-benefit: smooth observed-vs-predicted pairwise benefit.
#
# Returns Eavg and E50. Also returns the smoothed predictions for plotting.
# -----------------------------------------------------------------------------
compute_calibration_for_benefit <- function(tau_predicted, tau_observed) {
  
  fit <- stats::loess(tau_observed ~ tau_predicted, degree = 2)
  tau_smoothed <- fit$fitted
  
  Eavg <- mean(abs(tau_smoothed - tau_predicted))
  E50  <- stats::median(abs(tau_smoothed - tau_predicted))
  
  list(
    Eavg         = Eavg,
    E50          = E50,
    tau_smoothed = tau_smoothed
  )
}

# -----------------------------------------------------------------------------
# RMSE between learned predicted benefit and oracle ("th") predicted benefit.
# Both inputs are pairwise (averaged within matched pair).
# -----------------------------------------------------------------------------
compute_rmse_vs_oracle <- function(tau_predicted, tau_predicted_oracle) {
  sqrt(mean((tau_predicted_oracle - tau_predicted) ^ 2))
}

# -----------------------------------------------------------------------------
# Calibration-for-benefit plot.
# -----------------------------------------------------------------------------
plot_calibration_for_benefit <- function(tau_predicted, tau_smoothed,
                                         xlim = c(-0.1, 1),
                                         ylim = c(-0.1, 1)) {
  df <- data.frame(tau_predicted = tau_predicted,
                   tau_smoothed  = tau_smoothed)
  
  ggplot2::ggplot(df, ggplot2::aes(x = tau_predicted, y = tau_smoothed)) +
    ggplot2::geom_point(shape = 9, size = 2, color = "#feb24c") +
    ggplot2::geom_abline(intercept = 0, linetype = "dashed",
                         color = "#756bb1") +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(ylim[1], ylim[2]) +
    ggplot2::labs(x = "Predicted pairwise treatment effect",
                  y = "Observed pairwise treatment effect",
                  color = " ") +
    ggplot2::theme_classic()
}

save_calibration_plot <- function(plot_obj, dir_plot, filename = "LMPlot.pdf") {
  if (!dir.exists(dir_plot)) dir.create(dir_plot, recursive = TRUE)
  pdf(file.path(dir_plot, filename))
  print(plot_obj, newpage = FALSE)
  dev.off()
  invisible(NULL)
}