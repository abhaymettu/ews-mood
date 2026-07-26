## sim/02_null_calibration.R -- what does "significant" mean for a rolling-window
## Kendall tau?
##
## The applied EWS literature routinely reports a Kendall tau and its nominal
## p-value. That p-value assumes independent observations. Rolling windows with
## a 2-day step and a 30-day width share ~93% of their data with their
## neighbours, so the indicator series is enormously autocorrelated and the
## nominal test is anticonservative by a large factor.
##
## This script measures the factor, and produces the empirical thresholds that
## every later script uses. Nothing downstream is allowed to use a nominal
## p-value.
##
## It also checks the drift null does not accidentally tip -- if it does, we are
## conditioning on survival and the FPR is biased.

source("R/setup.R")
set.seed(SEED)

NREP  <- 400
SPAN  <- 180            # evaluated span, matching the longest cell of the power grid
WIN   <- c(20, 30, 45, 60)
INDS  <- c("variance", "ac1_naive", "ac1_ou")

sim_null <- function(i, scenario) {
  set.seed(SEED + 1000 * (scenario == "drift") + i)
  p <- simulate_latent(SIM_DAYS, scenario)
  ## only tips inside the evaluated span matter: a null that tips at day 250
  ## when we stop looking at day 190 is still a valid null for this test
  tt <- find_transition(p)
  tipped <- is.finite(tt) && tt <= REF_DAY * 24
  ## Observation is clean -- no measurement noise, no Likert, no missingness --
  ## because this script isolates the STATISTICS of the trend test, not realism.
  ## Sampling is at a realistic 10 prompts/day rather than continuously, so the
  ## number of windows and their overlap match the regime we will actually work
  ## in; that is what drives the null distribution.
  d <- sample_path(p, esm_times(SPAN, 10, day0 = (REF_DAY - SPAN) * 24))
  pre <- detrend(d$x, d$t, "gaussian", bw_days = 14)   # once: it is the slow step
  out <- expand.grid(win = WIN, indicator = INDS, stringsAsFactors = FALSE)
  out$tau <- mapply(function(w, k)
    ews_tau(NULL, d$t, k, win_days = w, step_days = 2, pre = pre),
    out$win, out$indicator)
  out$rep <- i; out$scenario <- scenario; out$tipped <- tipped
  out
}

ncore <- max(1L, parallel::detectCores() - 1L)
cat("running", NREP, "reps x 2 nulls on", ncore, "cores\n")
null_res <- do.call(rbind, c(
  parallel::mclapply(1:NREP, sim_null, scenario = "static", mc.cores = ncore),
  parallel::mclapply(1:NREP, sim_null, scenario = "drift",  mc.cores = ncore)
))

cat(sprintf("\ndrift nulls that tipped inside the evaluated span (day 0-%d): %.1f%%  (want 0 -- else we are
conditioning on survival and the false-alarm rate is biased)\n",
            REF_DAY, 100 * mean(subset(null_res, scenario == "drift" & win == 30 &
                                       indicator == "variance")$tipped)))

## --- empirical vs nominal ----------------------------------------------------
## nominal: two-sided Kendall tau test at alpha = .05 on the rolling series.
## We recover the number of windows to get the nominal critical value right.
nwin <- function(w, span_days = SPAN) floor((span_days - w) / 2) + 1
crit_nominal <- function(n) qnorm(0.95) * sqrt(2 * (2 * n + 5) / (9 * n * (n - 1)))

thr <- do.call(rbind, lapply(split(subset(null_res, scenario == "static"),
                                   ~ win + indicator), function(z) {
  n <- nwin(z$win[1])
  data.frame(win = z$win[1], indicator = z$indicator[1],
             thr_emp95 = quantile(z$tau, .95, na.rm = TRUE),
             thr_nominal = crit_nominal(n),
             n_windows = n,
             fpr_if_nominal = mean(z$tau > crit_nominal(n), na.rm = TRUE))
}))
rownames(thr) <- NULL
cat("\n--- static null: empirical vs nominal 5% critical value for tau ---\n")
print(thr, digits = 3)

## --- clinical null: how often does the calibrated test fire on a
##     non-tipping but genuinely slowing person? -------------------------------
drift_fpr <- do.call(rbind, lapply(split(subset(null_res, scenario == "drift"),
                                         ~ win + indicator), function(z) {
  tt <- thr$thr_emp95[thr$win == z$win[1] & thr$indicator == z$indicator[1]]
  data.frame(win = z$win[1], indicator = z$indicator[1],
             fpr_drift = mean(z$tau > tt, na.rm = TRUE))
}))
rownames(drift_fpr) <- NULL
cat("\n--- clinical (drift) null: alarm rate at the 5%-calibrated threshold ---\n")
print(merge(thr[, c("win", "indicator", "thr_emp95")], drift_fpr), digits = 3)

write.csv(thr, file.path(TAB, "02_thresholds.csv"), row.names = FALSE)
write.csv(null_res, file.path(TAB, "02_null_taus.csv"), row.names = FALSE)
saveRDS(list(thr = thr, drift_fpr = drift_fpr, null_res = null_res),
        file.path(TAB, "02_null_calibration.rds"))

p <- ggplot(null_res, aes(tau, fill = scenario)) +
  geom_density(alpha = .45, colour = NA) +
  geom_vline(data = thr, aes(xintercept = thr_nominal), linetype = 3, colour = "firebrick") +
  geom_vline(data = thr, aes(xintercept = thr_emp95), linetype = 2) +
  facet_grid(indicator ~ win, labeller = label_both) +
  scale_fill_manual(values = c(static = "grey40", drift = "darkorange")) +
  labs(title = "Null distribution of the rolling-window Kendall tau",
       subtitle = "dotted red = nominal 5% critical value; dashed black = empirical 5% from the static null",
       x = "Kendall tau", y = NULL)
ggsave(file.path(FIG, "02_null_distributions.png"), p, width = 10, height = 6, dpi = 150)
