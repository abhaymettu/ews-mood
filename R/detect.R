## R/detect.R -- early-warning-signal detectors
##
## Deliberately hand-rolled rather than taken from `earlywarnings`. Three
## reasons, in order of importance:
##   (a) `generic_ews` assumes an evenly-spaced ts object. ESM data is not
##       evenly spaced and pretending otherwise is exactly the artefact we are
##       here to quantify.
##   (b) the whole point of Phase 0 is to know what the detector does. A black
##       box you cannot defend in an interview is worse than 40 lines you can.
##   (c) every operation below is base R. cor(), optim(), lm().
## We still cross-check against `earlywarnings` on the ideal evenly-spaced case
## in sim/01_ideal.R -- if our AC1 disagrees with theirs on data they can handle,
## we are wrong.

## Detrending -----------------------------------------------------------------
## Non-negotiable step, and the single most contestable one (D4). Approaching
## the fold the *mean* slides from x*~1.15 toward x*~0.58. A rolling variance
## computed on undetrended data picks up that slide and reports "rising
## variance" even with a completely static noise process. So detrending is
## required to make the test about the residual dynamics rather than the trend.
##
## But detrend too hard and you remove the slow variance you were trying to
## measure. There is no correct bandwidth. Phase 2 sweeps it; Phase 0 shows how
## much it moves the answer.
##
## Gaussian kernel smoother (Nadaraya-Watson) works natively on irregular time,
## which is why it is the default here rather than a filter that assumes a grid.
detrend <- function(y, t, method = c("gaussian", "linear", "none"), bw_days = 14) {
  method <- match.arg(method)
  switch(method,
    none     = y - mean(y),
    linear   = residuals(lm(y ~ t)),
    gaussian = {
      ## Evaluated on a coarse grid and linearly interpolated back. The smoother
      ## has bandwidth of weeks, so it cannot bend appreciably between grid
      ## points spaced a few hours apart; this is a numerical shortcut, not an
      ## approximation to the estimator. Exact O(n^2) version is kept below and
      ## checked against this one in sim/01_ideal.R.
      bw <- bw_days * 24
      g  <- seq(min(t), max(t), length.out = min(length(t), 400L))
      fg <- vapply(g, function(gi) {
        w <- dnorm(t - gi, sd = bw); sum(w * y) / sum(w)
      }, numeric(1))
      y - approx(g, fg, xout = t, rule = 2)$y
    }
  )
}

## Indicator 1: variance -------------------------------------------------------
ind_var <- function(y, t) var(y)

## Indicator 2: naive lag-1 autocorrelation ------------------------------------
## "Naive" = treats consecutive *observations* as lag-1 regardless of the actual
## elapsed time between them. This is what almost every applied EWS paper does.
## On ESM data the gaps are a mixture of ~1.5h within-day and ~11h overnight, and
## after 25% missingness the mixture shifts. The estimate is a time-average of
## exp(lambda*gap) over a gap distribution that is itself moving. We keep it
## precisely so we can measure how badly it misleads.
ind_ac1_naive <- function(y, t) {
  if (length(y) < 10) return(NA_real_)
  cor(head(y, -1), tail(y, -1))
}

## Indicator 3: within-day AC1 --------------------------------------------------
## Cheap partial fix: use only consecutive pairs whose gap is close to the modal
## within-day spacing. Removes the overnight gap and the post-missingness long
## gaps. Costs sample size.
ind_ac1_withinday <- function(y, t, max_gap_h = 3) {
  d <- diff(t)
  ok <- d <= max_gap_h
  if (sum(ok) < 10) return(NA_real_)
  cor(head(y, -1)[ok], tail(y, -1)[ok])
}

## Indicator 4: continuous-time OU recovery rate -------------------------------
## The principled version. The linearised system near a stable point is an OU
## process, whose exact transition density over an arbitrary gap D is
##   y_{i+1} | y_i ~ N( mu + (y_i - mu) e^{-lambda D_i},  s2 (1 - e^{-2 lambda D_i}) )
## so irregular spacing is not a nuisance to be patched, it is just data. We fit
## (mu, lambda, s2) by exact ML over the observed pairs and report lambda per
## hour. CSD is lambda -> 0, so this estimates the theoretical quantity directly
## instead of a gap-contaminated proxy for it.
##
## Reported as an "AC1-equivalent" exp(-lambda * ref_gap) so it is on the same
## 0-1 scale as the other indicators and comparable across sampling densities.
##
## Known limitation: measurement noise biases lambda UPWARD (see below). This
## estimator does not separate the two; ind_ou_ar1ma1 does.
fit_ou <- function(y, t) {
  d <- diff(t)
  y0 <- head(y, -1); y1 <- tail(y, -1)
  ok <- is.finite(d) & d > 0
  if (sum(ok) < 15) return(c(lambda = NA_real_, s2 = NA_real_))
  d <- d[ok]; y0 <- y0[ok]; y1 <- y1[ok]

  nll <- function(p) {
    mu <- p[1]; lam <- exp(p[2]); s2 <- exp(p[3])
    e <- exp(-lam * d)
    v <- s2 * (1 - e^2)
    if (any(!is.finite(v)) || any(v <= 1e-12)) return(1e10)
    m <- mu + (y0 - mu) * e
    0.5 * sum(log(2 * pi * v) + (y1 - m)^2 / v)
  }
  st <- c(mean(y), log(0.3), log(max(var(y), 1e-4)))
  fit <- tryCatch(optim(st, nll, method = "BFGS",
                        control = list(maxit = 400, reltol = 1e-9)),
                  error = function(e) NULL)
  if (is.null(fit) || fit$convergence != 0) return(c(lambda = NA_real_, s2 = NA_real_))
  c(lambda = exp(fit$par[2]), s2 = exp(fit$par[3]))
}

ind_ou_ac1 <- function(y, t, ref_gap = 1.5) {
  p <- fit_ou(y, t)
  exp(-p["lambda"] * ref_gap)
}

## Indicator 5: OU + measurement error -----------------------------------------
## The econometric point. An AR(1) signal observed with additive white noise is
## an ARMA(1,1); regressing y_{t+1} on y_t gives an autocorrelation attenuated by
## the reliability ratio
##   plim(rho_hat) = rho * s2 / (s2 + sigma_m^2)
## As the system approaches the fold, s2 rises, so the attenuation *shrinks* --
## measurement error does not just add a constant bias, it adds a bias that
## moves in the same direction as the signal. Some of the "rising AC1" in a
## noisy ESM series is the reliability ratio improving, not the recovery rate
## slowing. The two are only separable if you model them separately, which is
## what this does: state-space OU with an explicit measurement-noise variance,
## fit by Kalman filter.
##
## Costs: 4 parameters on short windows, and it is the slowest indicator by an
## order of magnitude. Whether that cost buys anything is an empirical question
## and is answered in sim/03_power.R.
fit_ou_me <- function(y, t) {
  d <- c(NA, diff(t))
  n <- length(y)
  if (n < 25) return(c(lambda = NA_real_, s2 = NA_real_, sm2 = NA_real_))
  nll <- function(p) {
    mu <- p[1]; lam <- exp(p[2]); s2 <- exp(p[3]); sm2 <- exp(p[4])
    xh <- y[1] - mu; Ph <- s2                      # diffuse-ish init at stationarity
    ll <- 0
    for (i in 2:n) {
      e  <- exp(-lam * d[i])
      xp <- xh * e
      Pp <- Ph * e^2 + s2 * (1 - e^2)
      v  <- Pp + sm2
      if (!is.finite(v) || v <= 1e-12) return(1e10)
      r  <- (y[i] - mu) - xp
      ll <- ll + log(2 * pi * v) + r^2 / v
      K  <- Pp / v
      xh <- xp + K * r
      Ph <- Pp * (1 - K)
    }
    0.5 * ll
  }
  st <- c(mean(y), log(0.3), log(max(var(y) * 0.5, 1e-4)), log(max(var(y) * 0.5, 1e-4)))
  fit <- tryCatch(optim(st, nll, method = "Nelder-Mead",
                        control = list(maxit = 1500)),
                  error = function(e) NULL)
  if (is.null(fit) || fit$convergence != 0)
    return(c(lambda = NA_real_, s2 = NA_real_, sm2 = NA_real_))
  c(lambda = exp(fit$par[2]), s2 = exp(fit$par[3]), sm2 = exp(fit$par[4]))
}

ind_ou_me_ac1 <- function(y, t, ref_gap = 1.5) {
  exp(-fit_ou_me(y, t)["lambda"] * ref_gap)
}

INDICATORS <- list(
  variance     = ind_var,
  ac1_naive    = ind_ac1_naive,
  ac1_withinday = ind_ac1_withinday,
  ac1_ou       = ind_ou_ac1,
  ac1_ou_me    = ind_ou_me_ac1
)

## Rolling windows -------------------------------------------------------------
## Windows are defined in TIME, not in observation count. With missingness the
## two differ, and a count-based window silently stretches over more calendar
## time exactly when data is scarcest -- which is exactly when mood is worst
## under MNAR. Time-based windows keep the comparison honest.
rolling <- function(y, t, fun, win_days = 30, step_days = 2, min_n = 20, ...) {
  w  <- win_days * 24
  if (!length(t) || diff(range(t)) < w) return(NULL)   # window wider than the series
  st <- seq(min(t), max(t) - w, by = step_days * 24)
  if (!length(st)) return(NULL)
  out <- lapply(st, function(s) {
    ix <- t >= s & t < s + w
    if (sum(ix) < min_n) return(c(end = s + w, val = NA_real_, n = sum(ix)))
    c(end = s + w, val = as.numeric(fun(y[ix], t[ix], ...)), n = sum(ix))
  })
  as.data.frame(do.call(rbind, out))
}

## Trend statistic -------------------------------------------------------------
## Kendall's tau of the rolling indicator against window time. This is the field
## standard (Dakos et al. 2012).
##
## THE central methodological trap: the tau from overlapping rolling windows is
## massively autocorrelated, so its nominal p-value is meaningless -- it will
## declare significance on pure noise at rates far above 5%. Every threshold in
## this project is therefore calibrated against simulated null runs, never
## against the analytic null. sim/03_power.R does exactly that, and sim/02_null_calibration.R
## measures how wrong the nominal p-value is.
kendall_tau <- function(r) {
  if (is.null(r)) return(NA_real_)
  ok <- is.finite(r$val)
  if (sum(ok) < 5) return(NA_real_)
  suppressWarnings(cor(r$end[ok], r$val[ok], method = "kendall"))
}

## Exact (slow) reference implementation, used only to verify the grid version.
detrend_exact <- function(y, t, bw_days = 14) {
  bw <- bw_days * 24
  y - vapply(seq_along(t), function(i) {
    w <- dnorm(t - t[i], sd = bw); sum(w * y) / sum(w)
  }, numeric(1))
}

## Convenience: indicator -> rolling series -> tau, in one call.
## Pass `pre` (an already-detrended series) when looping over several indicators
## or window widths on the same data -- detrending is the expensive step and it
## does not depend on either.
ews_tau <- function(y, t, indicator, win_days = 30, step_days = 2,
                    detrend_method = "gaussian", bw_days = 14, pre = NULL) {
  yy <- if (is.null(pre)) detrend(y, t, detrend_method, bw_days) else pre
  kendall_tau(rolling(yy, t, INDICATORS[[indicator]], win_days, step_days))
}
