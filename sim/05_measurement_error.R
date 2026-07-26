## sim/05_measurement_error.R -- can a state-space estimator separate true
## slowing from shrinking attenuation?
##
## The problem, from sim/03: an AR(1) observed with additive white noise is an
## ARMA(1,1), and the naive lag-1 autocorrelation is attenuated by the
## reliability ratio s2/(s2 + sigma_m^2). Since s2 rises as the fold is
## approached, the attenuation shrinks along the approach and the naive AC1
## rises for a reason unrelated to the recovery rate. Measured there: ~59% of
## the apparent rise at sigma_m = 0.4.
##
## That is a bias pointing the same way as the hypothesis, in a literature that
## tests the hypothesis by looking for a rise. So it is worth asking whether the
## fix works.
##
## `ind_ou_me_ac1` fits a continuous-time OU state process with an explicit
## measurement-noise variance by Kalman filter, so lambda is identified
## separately from sigma_m. Two questions:
##   Q1  does it recover sigma_m and lambda when we know the truth?
##   Q2  does it beat the moment-based indicators on discrimination (AUC against
##       the clinical null), and is that worth ~50x the compute?
##
## Missingness is off and the design is fixed at 180 days / 10 prompts so that
## measurement error is the only thing varying.

source("R/setup.R")

NREP     <- 150
DAYS     <- 180
PPD      <- 10
SIGMA_M  <- c(0.1, 0.4, 0.8)
INDS     <- c("ac1_naive", "ac1_ou", "ac1_ou_me", "variance")

## --- Q1: parameter recovery --------------------------------------------------
## Simulate a stationary OU with known lambda and known measurement noise on the
## actual ESM sampling grid, then see what the estimators return. If the
## state-space fit cannot recover the truth here it cannot be trusted on the
## curved case.
cat("--- Q1: parameter recovery on a stationary OU at the ESM sampling grid ---\n")
recov <- do.call(rbind, parallel::mclapply(1:200, function(i) {
  set.seed(SEED + 300000 + i)
  lam <- 0.6                                     # per hour, ~ the euthymic branch
  s2  <- 0.05
  tt  <- esm_times(45, PPD)
  d   <- diff(tt)
  x   <- numeric(length(tt)); x[1] <- rnorm(1, 0, sqrt(s2))
  for (j in 2:length(tt)) {
    e <- exp(-lam * d[j - 1])
    x[j] <- x[j - 1] * e + rnorm(1, 0, sqrt(s2 * (1 - e^2)))
  }
  do.call(rbind, lapply(SIGMA_M, function(sm) {
    y  <- x + rnorm(length(x), 0, sm)
    p1 <- fit_ou(y, tt)
    p2 <- fit_ou_me(y, tt)
    data.frame(sigma_m = sm,
               lam_ou    = unname(p1["lambda"]),
               lam_ou_me = unname(p2["lambda"]),
               sm_hat    = sqrt(unname(p2["sm2"])),
               ac1_naive = cor(head(y, -1), tail(y, -1)))
  }))
}, mc.cores = max(1L, parallel::detectCores() - 1L)))

rc <- aggregate(cbind(lam_ou, lam_ou_me, sm_hat, ac1_naive) ~ sigma_m, recov,
                function(z) mean(z, na.rm = TRUE))
rc$lam_true <- 0.6
rc$ac1_true <- exp(-0.6 * mean(diff(esm_times(45, PPD))))
cat("true lambda = 0.60 per hour\n")
print(rc, row.names = FALSE, digits = 3)
cat("\n-> ac1_naive collapses toward 0 as sigma_m grows (attenuation);",
    "\n   fit_ou inherits the same bias in lambda (biased UPWARD: it reads",
    "\n   measurement noise as fast mean reversion);",
    "\n   fit_ou_me should recover both lambda and sigma_m.\n")
write.csv(rc, file.path(TAB, "05_recovery.csv"), row.names = FALSE)

## --- Q2: discrimination ------------------------------------------------------
run_rep <- function(i, scenario) {
  set.seed(SEED + 100000 * match(scenario, c("transition", "static", "drift")) + i)
  p  <- simulate_latent(SIM_DAYS, scenario)
  tt <- find_transition(p)
  if (scenario == "transition") {
    if (is.na(tt) || tt < DAYS * 24) return(NULL)
    end_h <- tt
  } else {
    if (!is.na(tt) && tt <= REF_DAY * 24) return(NULL)
    end_h <- REF_DAY * 24
  }
  do.call(rbind, lapply(SIGMA_M, function(sm) {
    set.seed(SEED + 900000 + i)
    o <- observe(p, end_h, DAYS, PPD, sigma_m = sm, K = 7)
    if (is.null(o)) return(NULL)
    pre <- detrend(o$y, o$t, "gaussian", bw_days = DAYS / 8)
    data.frame(rep = i, scenario = scenario, sigma_m = sm, indicator = INDS,
               tau = vapply(INDS, function(k)
                 ews_tau(NULL, o$t, k, win_days = DAYS / 2, step_days = 6, pre = pre),
                 numeric(1)))
  }))
}

cat("\n--- Q2: discrimination (this is the slow part) ---\n")
t0 <- Sys.time()
res <- do.call(rbind, unlist(lapply(c("transition", "static", "drift"), function(s)
  parallel::mclapply(1:NREP, run_rep, scenario = s,
                     mc.cores = max(1L, parallel::detectCores() - 1L))),
  recursive = FALSE))
cat(sprintf("done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

perf <- do.call(rbind, lapply(split(res, res[c("sigma_m", "indicator")], drop = TRUE), function(z) {
  tr <- na.omit(z$tau[z$scenario == "transition"])
  st <- na.omit(z$tau[z$scenario == "static"])
  dr <- na.omit(z$tau[z$scenario == "drift"])
  if (length(st) < 30 || length(tr) < 30) return(NULL)
  thr <- quantile(st, 0.95)
  data.frame(sigma_m = z$sigma_m[1], indicator = z$indicator[1],
             thr = as.numeric(thr), power = mean(tr > thr),
             fpr_drift = mean(dr > thr),
             auc = mean(outer(tr, dr, ">")) + 0.5 * mean(outer(tr, dr, "==")),
             na_rate = mean(is.na(z$tau[z$scenario == "transition"])))
}))
rownames(perf) <- NULL
perf <- perf[order(perf$sigma_m, -perf$auc), ]
cat("\n--- power, false alarms and AUC by measurement-noise level ---\n")
print(perf, row.names = FALSE, digits = 3)
write.csv(perf, file.path(TAB, "05_measurement_error.csv"), row.names = FALSE)
saveRDS(list(recov = recov, perf = perf, res = res), file.path(TAB, "05_me.rds"))

p <- ggplot(perf, aes(factor(sigma_m), auc, fill = indicator)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(title = "Does modelling measurement error explicitly help?",
       subtitle = "AUC for discriminating a transition from a non-tipping deterioration; 0.5 = useless",
       x = "measurement noise (Likert points)", y = "AUC", fill = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "05_measurement_error.png"), p, width = 8, height = 5, dpi = 150)
