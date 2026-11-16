# =============================================================================
# src/metrics/decision_accuracy.R
# Population Benefit (Efthimiou et al. 2025, Stat Med).
#
# Groups (Figure 1):
#   G1: pred_ben > 0 & trt == 1     (model says treat, treated)
#   G2: pred_ben > 0 & trt == 0     (model says treat, untreated)
#   G3: pred_ben < 0 & trt == 1     (model says don't treat, treated)
#   G4: pred_ben < 0 & trt == 0     (model says don't treat, untreated)
#   agreed = 1  <=> G1 ∪ G4   (follow model)
#   agreed = 0  <=> G2 ∪ G3   (against model)
# =============================================================================
pb_survprob <- function(pred_ben, trt01, time, status, t.star) {
  
  N      <- length(pred_ben)
  agreed <- as.numeric((trt01 == 1) == (pred_ben > 0))
  
  d <- data.frame(time = time, status = status, trt01 = trt01,
                  pred_ben = pred_ben, agreed = agreed)
  
  # KM fits for the four unions ----------------------------------------------
  km14 <- survival::survfit(survival::Surv(time, status) ~ 1, data = d[d$agreed == 1, ])
  km23 <- survival::survfit(survival::Surv(time, status) ~ 1, data = d[d$agreed == 0, ])
  km13 <- survival::survfit(survival::Surv(time, status) ~ 1, data = d[d$trt01  == 1, ])
  km24 <- survival::survfit(survival::Surv(time, status) ~ 1, data = d[d$trt01  == 0, ])
  
  # KM fits for the four individual groups -----------------------------------
  km1 <- survival::survfit(survival::Surv(time, status) ~ 1,
                           data = d[d$trt01 == 1 & d$pred_ben > 0, ])
  km2 <- survival::survfit(survival::Surv(time, status) ~ 1,
                           data = d[d$trt01 == 0 & d$pred_ben > 0, ])
  km3 <- survival::survfit(survival::Surv(time, status) ~ 1,
                           data = d[d$trt01 == 1 & d$pred_ben < 0, ])
  km4 <- survival::survfit(survival::Surv(time, status) ~ 1,
                           data = d[d$trt01 == 0 & d$pred_ben < 0, ])
  
  # S(t*) and SE for each fit ------------------------------------------------
  s14 <- summary(km14, times = t.star, extend = TRUE)
  s23 <- summary(km23, times = t.star, extend = TRUE)
  s13 <- summary(km13, times = t.star, extend = TRUE)
  s24 <- summary(km24, times = t.star, extend = TRUE)
  s1  <- summary(km1,  times = t.star, extend = TRUE)
  s2  <- summary(km2,  times = t.star, extend = TRUE)
  s3  <- summary(km3,  times = t.star, extend = TRUE)
  s4  <- summary(km4,  times = t.star, extend = TRUE)
  
  # ---- PB : follow model vs against model ----------------------------------
  PB    <- s14$surv[1] - s23$surv[1]
  SE.PB <- sqrt(s14$std.err[1]^2 + s23$std.err[1]^2)
  
  # ---- PB_treat_all : ATE on survival probability scale --------------------
  PB.all    <- s13$surv[1] - s24$surv[1]
  SE.PB.all <- sqrt(s13$std.err[1]^2 + s24$std.err[1]^2)
  
  # ---- PB_MvT : follow model vs treat everyone -----------------------------
  p.neg   <- mean(pred_ben < 0)
  v.p.neg <- p.neg * (1 - p.neg) / N
  d43     <- s4$surv[1] - s3$surv[1]
  v43     <- s4$std.err[1]^2 + s3$std.err[1]^2
  
  PB.MvT    <- p.neg * d43
  SE.PB.MvT <- sqrt(v.p.neg * v43 + (p.neg^2) * v43 + (d43^2) * v.p.neg)
  
  # ---- PB_MvC : follow model vs treat no one -------------------------------
  p.pos   <- 1 - p.neg
  v.p.pos <- v.p.neg
  d12     <- s1$surv[1] - s2$surv[1]
  v12     <- s1$std.err[1]^2 + s2$std.err[1]^2
  
  PB.MvC    <- p.pos * d12
  SE.PB.MvC <- sqrt(v.p.pos * v12 + (p.pos^2) * v12 + (d12^2) * v.p.pos)
  
  # ---- pack ----------------------------------------------------------------
  c(
    PB.est           = PB,        PB.se           = SE.PB,
    PB.lo            = PB     - 1.96 * SE.PB,
    PB.hi            = PB     + 1.96 * SE.PB,
    
    PB_MvT.est       = PB.MvT,    PB_MvT.se       = SE.PB.MvT,
    PB_MvT.lo        = PB.MvT - 1.96 * SE.PB.MvT,
    PB_MvT.hi        = PB.MvT + 1.96 * SE.PB.MvT,
    
    PB_MvC.est       = PB.MvC,    PB_MvC.se       = SE.PB.MvC,
    PB_MvC.lo        = PB.MvC - 1.96 * SE.PB.MvC,
    PB_MvC.hi        = PB.MvC + 1.96 * SE.PB.MvC,
    
    PB_treat_all.est = PB.all,    PB_treat_all.se = SE.PB.all,
    PB_treat_all.lo  = PB.all - 1.96 * SE.PB.all,
    PB_treat_all.hi  = PB.all + 1.96 * SE.PB.all
  )
}

# -----------------------------------------------------------------------------
# PB family on the RMST scale.
# -----------------------------------------------------------------------------
pb_rmst <- function(pred_ben, trt01, time, status, t.star) {
  
  N      <- length(pred_ben)
  agreed <- as.numeric((trt01 == 1) == (pred_ben > 0))
  
  d <- data.frame(time = time, status = status, trt01 = trt01,
                  pred_ben = pred_ben, agreed = agreed)
  
  # Returns NA if the subset is too small or t.star is below smallest time.
  rmst_of <- function(idx) {
    if (sum(idx) < 2 || sum(d$status[idx]) < 1)  return(c(NA_real_, NA_real_))
    if (t.star < min(d$time[idx]))               return(c(NA_real_, NA_real_))
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = d[idx, ])
    tab <- tryCatch(summary(fit, rmean = t.star)$table, error = function(e) NULL)
    if (is.null(tab))                            return(c(NA_real_, NA_real_))
    rmean_name <- grep("rmean$",     names(tab), value = TRUE)[1]
    se_name    <- grep("se\\(rmean", names(tab), value = TRUE)[1]
    if (is.na(rmean_name) || is.na(se_name))     return(c(NA_real_, NA_real_))
    c(unname(tab[rmean_name]), unname(tab[se_name]))
  }
  
  r14 <- rmst_of(d$agreed == 1)
  r23 <- rmst_of(d$agreed == 0)
  r13 <- rmst_of(d$trt01  == 1)
  r24 <- rmst_of(d$trt01  == 0)
  r1  <- rmst_of(d$trt01 == 1 & d$pred_ben > 0)
  r2  <- rmst_of(d$trt01 == 0 & d$pred_ben > 0)
  r3  <- rmst_of(d$trt01 == 1 & d$pred_ben < 0)
  r4  <- rmst_of(d$trt01 == 0 & d$pred_ben < 0)
  
  # ---- PB : follow model vs against model ----------------------------------
  PB    <- r14[1] - r23[1]
  SE.PB <- sqrt(r14[2]^2 + r23[2]^2)
  
  # ---- PB_treat_all : ATE on RMST scale ------------------------------------
  PB.all    <- r13[1] - r24[1]
  SE.PB.all <- sqrt(r13[2]^2 + r24[2]^2)
  
  # ---- PB_MvT : follow model vs treat everyone -----------------------------
  p.neg   <- mean(pred_ben < 0)
  v.p.neg <- p.neg * (1 - p.neg) / N
  d43     <- r4[1] - r3[1]
  v43     <- r4[2]^2 + r3[2]^2
  PB.MvT    <- p.neg * d43
  SE.PB.MvT <- sqrt(v.p.neg * v43 + (p.neg^2) * v43 + (d43^2) * v.p.neg)
  
  # ---- PB_MvC : follow model vs treat no one -------------------------------
  p.pos   <- 1 - p.neg
  v.p.pos <- v.p.neg
  d12     <- r1[1] - r2[1]
  v12     <- r1[2]^2 + r2[2]^2
  PB.MvC    <- p.pos * d12
  SE.PB.MvC <- sqrt(v.p.pos * v12 + (p.pos^2) * v12 + (d12^2) * v.p.pos)
  
  # ---- pack ----------------------------------------------------------------
  c(
    PB_RMST.est           = PB,        PB_RMST.se           = SE.PB,
    PB_RMST.lo            = PB     - 1.96 * SE.PB,
    PB_RMST.hi            = PB     + 1.96 * SE.PB,
    
    PB_RMST_MvT.est       = PB.MvT,    PB_RMST_MvT.se       = SE.PB.MvT,
    PB_RMST_MvT.lo        = PB.MvT - 1.96 * SE.PB.MvT,
    PB_RMST_MvT.hi        = PB.MvT + 1.96 * SE.PB.MvT,
    
    PB_RMST_MvC.est       = PB.MvC,    PB_RMST_MvC.se       = SE.PB.MvC,
    PB_RMST_MvC.lo        = PB.MvC - 1.96 * SE.PB.MvC,
    PB_RMST_MvC.hi        = PB.MvC + 1.96 * SE.PB.MvC,
    
    PB_RMST_treat_all.est = PB.all,    PB_RMST_treat_all.se = SE.PB.all,
    PB_RMST_treat_all.lo  = PB.all - 1.96 * SE.PB.all,
    PB_RMST_treat_all.hi  = PB.all + 1.96 * SE.PB.all
  )
}

# -----------------------------------------------------------------------------
# Both scales together
# -----------------------------------------------------------------------------
pb_all <- function(pred_ben, trt01, time, status, t.star) {
  c(pb_survprob(pred_ben, trt01, time, status, t.star),
    pb_rmst    (pred_ben, trt01, time, status, t.star))
}