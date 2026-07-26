## R/model.R -- stochastic double-well model of mood with a drifting control parameter
##
## Model
## -----
##   dx = (-x^3 + x + c(t)) dt + sigma_p dW
##
## This is the gradient of the potential V(x) = x^4/4 - x^2/2 - c*x. It is the
## normal form of the cusp catastrophe and the standard toy model behind the
## critical-slowing-down literature (Scheffer et al. 2009; van de Leemput et al.
## PNAS 2014). x is a *latent* mood valence: x > 0 is read as the euthymic
## attractor, x < 0 as the depressed attractor. c is the slow control parameter
## (cumulative stress, medication dose, whatever you like).
##
## Fixed points solve x^3 - x - c = 0. Two stable branches exist for
## |c| < c_crit = 2/(3*sqrt(3)) = 0.3849. Lowering c past -c_crit annihilates the
## euthymic branch in a fold (saddle-node) bifurcation and the system drops into
## the depressed well. That fold is the "critical transition".
##
## Why CSD is expected: the local recovery rate at a stable point x* is
##   lambda(x*) = dF/dx = -3*x*^2 + 1
## On the euthymic branch, x* -> 1/sqrt(3) as c -> -c_crit, so lambda -> 0. The
## linearised process is Ornstein-Uhlenbeck with
##   AR(1) over gap D = exp(lambda * D)          -> 1
##   stationary variance = sigma_p^2 / (2|lambda|) -> Inf
## Rising autocorrelation and rising variance are therefore *derived*, not
## assumed. That derivation is also the yardstick: our detectors have to
## recover it on clean data or they are broken.
##
## CONTESTABLE (see DECISIONS.md D1, D2, D3): the cubic form, the linear ramp in
## c, and the timescale calibration are all choices, not facts about mood.

C_CRIT <- 2 / (3 * sqrt(3))

## Time calibration -----------------------------------------------------------
## The model is dimensionless. One model time unit (tu) has to be pinned to real
## hours or nothing else in the study means anything.
##
## Anchor: far from the bifurcation (c ~ 0.35, x* ~ 1.15) the recovery rate is
## lambda ~ -3*1.15^2 + 1 = -2.97 per tu. Empirical ESM affect series show lag-1
## autocorrelation ~0.4 at ~90-minute spacing ("emotional inertia"; Kuppens et
## al. 2010). Solving exp(-2.97 * 1.5/HPT) = 0.4 gives HPT ~ 4.9 hours per tu.
## We round to 5. CONTESTABLE (D3): 0.4 is a central value from a wide literature
## and the anchor scales every result about series length.
HOURS_PER_TU <- 5

## Analytic recovery rate on the upper (euthymic) branch, per hour.
## Used as ground truth to check the estimators against.
lambda_true <- function(cval, hours_per_tu = HOURS_PER_TU) {
  vapply(cval, function(cc) {
    if (cc <= -C_CRIT) return(NA_real_)          # branch no longer exists
    r <- Re(polyroot(c(-cc, -1, 0, 1)))          # roots of x^3 - x - c
    r <- r[abs(Im(polyroot(c(-cc, -1, 0, 1)))) < 1e-8]
    xs <- max(r)                                 # upper stable branch
    (-3 * xs^2 + 1) / hours_per_tu
  }, numeric(1))
}

## Simulate one latent path ---------------------------------------------------
## Euler-Maruyama. dt is in model time units; see sim/01_ideal.R for the
## step-size convergence check that justifies the default.
##
## scenario: we need TWO nulls, because they answer different questions (D7).
##   "transition" - c ramps linearly from c_start to c_end (< -c_crit): tips.
##   "static"     - c held at c_start. No change in the system at all. This is
##                  the CALIBRATION null: it isolates the purely statistical
##                  false-positive rate of a Kendall tau computed on overlapping
##                  rolling windows. All decision thresholds are set here.
##   "drift"      - c descends over the same span but ends `drift_gap` short of
##                  the fold. The system genuinely slows down and never tips.
##                  This is the CLINICAL null: the person who deteriorates,
##                  gets flagged, and does not become ill. Alarms here are false
##                  alarms in every sense that matters to a clinician even though
##                  real CSD is present. Reporting only the static-null FPR
##                  flatters the method badly; we report both.
##                  drift_gap = 0.35 was calibrated empirically (0/120 runs tip
##                  in 300 days). It cannot be pushed much lower: see D7.
## hold_days: c is held at c_start for this long before the ramp begins. Its
## only purpose is to push the transition late enough that every series length
## we want to study has that much history available in front of it. Without it,
## a 180-day window can only be cut from paths that happened to tip after day
## 180 -- and those are exactly the paths where noise did NOT tip the system
## early, i.e. the ones that got closest to the fold and slowed down most. That
## would make "longer series detect better" partly an artefact of conditioning
## on the easiest cases. See DECISIONS.md D13.
simulate_latent <- function(days,
                            scenario   = c("transition", "drift", "static"),
                            sigma_p    = 0.25,
                            c_start    = 0.35,
                            c_end      = -0.50,
                            drift_gap  = 0.35,   # null stops this far above the fold (calibrated: 0% tip)
                            hold_days  = if (exists("HOLD_DAYS")) HOLD_DAYS else 0,
                            dt         = 0.01,
                            hours_per_tu = HOURS_PER_TU,
                            x0         = NULL) {
  scenario <- match.arg(scenario)
  total_tu <- days * 24 / hours_per_tu
  n <- as.integer(ceiling(total_tu / dt))
  nh <- min(n, as.integer(hold_days * 24 / hours_per_tu / dt))   # held samples
  nr <- n - nh                                                   # ramping samples

  ramp <- function(to) c(rep(c_start, nh), seq(c_start, to, length.out = nr))
  cpath <- switch(scenario,
    transition = ramp(c_end),
    ## descends monotonically, ending drift_gap short of the fold. NOT capped
    ## and held there: a system parked near the fold escapes by noise with
    ## near-certainty over a study-length horizon (Kramers rate at c = -0.235 is
    ## 0.15 per tu, i.e. minutes-to-days), so "sits just short of tipping for
    ## months" is not a state that exists. See DECISIONS.md D7.
    drift      = ramp(-C_CRIT + drift_gap),
    static     = rep(c_start, n)
  )

  if (is.null(x0)) {
    ## start on the upper branch at c_start, not at an arbitrary value, so we
    ## are not measuring a transient relaxation as if it were dynamics
    r <- polyroot(c(-cpath[1], -1, 0, 1))
    x0 <- max(Re(r[abs(Im(r)) < 1e-8]))
  }

  x <- numeric(n)
  x[1] <- x0
  sq <- sigma_p * sqrt(dt)
  z <- rnorm(n, 0, sq)
  for (i in 2:n) {
    xi <- x[i - 1]
    x[i] <- xi + (-xi^3 + xi + cpath[i - 1]) * dt + z[i]
  }

  list(
    t_hours = (seq_len(n) - 1) * dt * hours_per_tu,
    x       = x,
    c       = cpath,
    dt      = dt,
    hours_per_tu = hours_per_tu,
    scenario = scenario
  )
}

## Locate the transition ------------------------------------------------------
## "First time the system is in the lower well and stays there." Requiring
## persistence stops us from calling a brief noise excursion a transition, which
## would poison every lead-time estimate.
find_transition <- function(path, thresh = -0.5, sustain_hours = 48) {
  below <- path$x < thresh
  if (!any(below)) return(NA_real_)
  step_h <- path$dt * path$hours_per_tu
  k <- max(1L, as.integer(sustain_hours / step_h))
  ## run-length trick: first index starting a run of >= k TRUEs
  cs <- c(0, cumsum(below))
  ok <- which((cs[(k + 1):length(cs)] - cs[1:(length(cs) - k)]) == k)
  if (!length(ok)) return(NA_real_)
  path$t_hours[ok[1]]
}
