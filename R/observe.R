## R/observe.R -- turn a continuous latent path into something that looks like ESM data
##
## Everything the real world does to the signal happens here, and each step is a
## separate switch so the power analysis can attribute damage to a specific
## cause rather than to "realism" in general.
##
##   1. discrete prompts inside a waking window, with an overnight gap
##   2. measurement noise (independent of process noise -- this distinction
##      turns out to matter more than anything else, see DECISIONS.md D5)
##   3. Likert discretisation with floor/ceiling
##   4. missingness, MCAR or mood-dependent (MNAR)

## Prompt schedule ------------------------------------------------------------
## Block-randomised, which is what real ESM protocols do: split the waking
## window into equal blocks and place one prompt at random inside the middle
## `jitter` fraction of each block. This guarantees a minimum spacing (people
## cannot be beeped twice in five minutes) while still being unpredictable.
esm_times <- function(n_days, prompts_per_day, wake = c(9, 22), jitter = 0.6,
                      day0 = 0) {
  span  <- diff(wake) / prompts_per_day
  edges <- wake[1] + (0:(prompts_per_day - 1)) * span
  lo    <- edges + span * (1 - jitter) / 2
  out   <- vapply(seq_len(n_days) - 1L, function(d) {
    d * 24 + lo + runif(prompts_per_day, 0, span * jitter)
  }, numeric(prompts_per_day))
  sort(as.vector(out)) + day0
}

## Sample the latent path at the prompt times (nearest fine-grid point;
## the grid step is ~3 min so interpolation error is negligible).
sample_path <- function(path, times) {
  step_h <- path$dt * path$hours_per_tu
  idx <- pmin(pmax(round(times / step_h) + 1L, 1L), length(path$x))
  data.frame(t = times, x = path$x[idx], c = path$c[idx])
}

## Response model -------------------------------------------------------------
## latent -> reported item score.
##   scale_pts : Likert points per unit of latent x. With the default well
##               separation of ~2.3 latent units and scale_pts = 2, the euthymic
##               state sits near 6.3 and the depressed state near 1.7 on a 1-7
##               item: a full swing that uses the scale without pinning it.
##   sigma_m   : measurement / momentary-situational noise, in Likert points,
##               added BEFORE rounding. This is white: it does not autocorrelate
##               and it does not grow near the bifurcation. It is pure dilution.
##   K         : number of response categories (0 = leave continuous, for
##               isolating the effect of discretisation).
respond <- function(x, sigma_m = 0.3, scale_pts = 2, center = 4, K = 7) {
  y <- center + scale_pts * x + rnorm(length(x), 0, sigma_m)
  if (K > 0) y <- pmin(pmax(round(y), 1), K)
  y
}

## Missingness ----------------------------------------------------------------
## MNAR mechanism: P(miss) rises as mood falls, because people skip prompts when
## they feel worst. Both mechanisms are calibrated to the SAME marginal missing
## rate, so any difference between them is attributable to the mechanism and not
## to sample size. (D6)
apply_missing <- function(obs, mechanism = c("none", "mcar", "mnar"),
                          rate = 0.25, slope = 1.5) {
  mechanism <- match.arg(mechanism)
  n <- nrow(obs)
  if (mechanism == "none" || rate <= 0) return(obs)
  keep <- if (mechanism == "mcar") {
    runif(n) > rate
  } else {
    ## solve intercept so mean(plogis(a - slope*x)) == rate
    f <- function(a) mean(plogis(a - slope * obs$x)) - rate
    a <- tryCatch(uniroot(f, c(-20, 20))$root, error = function(e) qlogis(rate))
    runif(n) > plogis(a - slope * obs$x)
  }
  obs[keep, , drop = FALSE]
}

## Full pipeline --------------------------------------------------------------
## Returns the observed series for one configuration. `window_days` takes the
## last N days ending at `end_hour` -- for transition runs we set end_hour to the
## transition time, so "180 days of data" means the 180 days immediately before
## the event, which is the situation a clinician is actually in.
observe <- function(path, end_hour, window_days, prompts_per_day,
                    sigma_m = 0.3, K = 7, scale_pts = 2,
                    mechanism = "none", miss_rate = 0.25,
                    wake = c(9, 22)) {
  start_hour <- end_hour - window_days * 24
  if (start_hour < 0) return(NULL)
  ## align the schedule to whole days so the overnight gap is where it should be
  day0 <- floor(start_hour / 24) * 24
  tt <- esm_times(window_days + 1, prompts_per_day, wake = wake, day0 = day0)
  tt <- tt[tt >= start_hour & tt <= end_hour]
  if (length(tt) < 20) return(NULL)

  obs <- sample_path(path, tt)
  obs$y <- respond(obs$x, sigma_m = sigma_m, scale_pts = scale_pts, K = K)
  obs <- apply_missing(obs, mechanism, miss_rate)
  obs$t_day <- obs$t / 24
  obs
}
