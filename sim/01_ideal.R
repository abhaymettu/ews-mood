## sim/01_ideal.R -- Phase 0 step 2: does CSD appear on clean, dense, evenly
## spaced data, where theory says it must?
##
## This is a gate, not a result. If AC1 and variance do not rise here, the model
## or the detector is broken and nothing downstream means anything.
##
## Three checks:
##   A. Euler-Maruyama step-size convergence (is dt = 0.01 fine enough?)
##   B. does the empirical rolling AC1/variance track the analytic lambda(c)?
##   C. sanity: no rise in the non-tipping null.

source("R/setup.R")
set.seed(SEED)

DAYS <- SIM_DAYS

## --- A. step-size convergence ------------------------------------------------
## Same noise realisation cannot be shared across dt, so compare distributions:
## 40 paths at each dt, compare the transition-time distribution and the
## stationary variance early in the run. If dt = 0.01 matches dt = 0.002 we use
## the cheap one.
conv <- do.call(rbind, lapply(c(0.02, 0.01, 0.002), function(dd) {
  set.seed(SEED + 1)
  r <- t(vapply(1:40, function(i) {
    p <- simulate_latent(DAYS, "transition", dt = dd)
    early <- p$x[p$t_hours < 24 * (HOLD_DAYS - 10)]
    c(tt = find_transition(p) / 24, sd_early = sd(early), mean_early = mean(early))
  }, numeric(3)))
  data.frame(dt = dd, tt = mean(r[, "tt"], na.rm = TRUE),
             tt_sd = sd(r[, "tt"], na.rm = TRUE),
             sd_early = mean(r[, "sd_early"]), mean_early = mean(r[, "mean_early"]))
}))
cat("\n--- A. Euler-Maruyama step-size convergence ---\n"); print(conv, digits = 4)
write.csv(conv, file.path(TAB, "01_stepsize_convergence.csv"), row.names = FALSE)

## --- B. the gate: CSD on clean data ------------------------------------------
## "Clean" = dense (every 30 min), no measurement noise, no missingness, no
## Likert rounding. The only realism is that we observe discretely.
set.seed(SEED + 7)
path <- simulate_latent(DAYS, "transition")
t_tr <- find_transition(path)
stopifnot("no transition occurred; check c_end / sigma_p" = is.finite(t_tr))
cat(sprintf("\ntransition at day %.1f (c = %.3f, c_crit = %.3f)\n",
            t_tr / 24, path$c[round(t_tr / (path$dt * path$hours_per_tu))], C_CRIT))

grid_t <- seq(0, t_tr, by = 0.5)                    # every 30 min
d <- sample_path(path, grid_t)
d$y <- d$x                                          # no measurement model yet

res <- lapply(c("variance", "ac1_naive", "ac1_ou"), function(k) {
  r <- rolling(detrend(d$y, d$t, "gaussian", bw_days = 14), d$t, INDICATORS[[k]],
               win_days = 40, step_days = 2)
  r$indicator <- k; r
})
res <- do.call(rbind, res)
res$day <- res$end / 24

taus <- vapply(split(res, res$indicator), function(z) kendall_tau(z), numeric(1))
cat("\n--- B. Kendall tau on clean data (want all clearly positive) ---\n")
print(round(taus, 3))

## analytic prediction: AC1-equivalent implied by lambda(c) at each window end
res$lambda_true <- NA_real_
step_h <- path$dt * path$hours_per_tu
cc <- path$c[pmin(round(res$end / step_h) + 1, length(path$c))]
res$pred_ac1 <- exp(-abs(lambda_true(cc)) * 1.5)

## --- C. null: same machinery, no approaching fold ----------------------------
set.seed(SEED + 11)
pnull <- simulate_latent(DAYS, "drift")
cat(sprintf("null path: c ranges %.3f to %.3f (fold at %.3f), transition = %s\n",
            min(pnull$c), max(pnull$c), -C_CRIT,
            ifelse(is.na(find_transition(pnull)), "none (good)", "TIPPED (bad)")))
dn <- sample_path(pnull, seq(0, t_tr, by = 0.5)); dn$y <- dn$x
taus_null <- vapply(c("variance", "ac1_naive", "ac1_ou"), function(k)
  ews_tau(dn$y, dn$t, k, win_days = 40, step_days = 2), numeric(1))
cat("\n--- C. Kendall tau on the non-tipping null (single run) ---\n")
print(round(taus_null, 3))

## --- figure ------------------------------------------------------------------
lab <- c(variance = "Variance", ac1_naive = "AC1 (naive lag-1)",
         ac1_ou = "AC1-equiv (continuous-time OU)")
res$indicator <- factor(res$indicator, names(lab), lab)

p_state <- ggplot(data.frame(day = path$t_hours / 24, x = path$x, c = path$c)) +
  geom_line(aes(day, x), linewidth = .2, colour = "grey30") +
  geom_vline(xintercept = t_tr / 24, colour = "firebrick", linetype = 2) +
  labs(title = "Latent mood, stochastic double well with a slowly falling control parameter",
       subtitle = sprintf("transition at day %.0f; dashed line = tipping point", t_tr / 24),
       x = NULL, y = "latent mood x")

p_ind <- ggplot(res, aes(day, val)) +
  geom_line(colour = "steelblue4") +
  geom_line(data = subset(res, indicator != "Variance"),
            aes(day, pred_ac1), colour = "firebrick", linetype = 3) +
  geom_vline(xintercept = t_tr / 24, colour = "firebrick", linetype = 2) +
  facet_wrap(~indicator, scales = "free_y", ncol = 1) +
  labs(x = "day", y = NULL,
       caption = "40-day rolling windows, Gaussian detrend (bw 14d). Dotted red = analytic prediction from lambda(c).")

ggsave(file.path(FIG, "01_ideal_state.png"), p_state, width = 8, height = 3, dpi = 150)
ggsave(file.path(FIG, "01_ideal_indicators.png"), p_ind, width = 8, height = 7, dpi = 150)

saveRDS(list(taus = taus, taus_null = taus_null, res = res, t_tr = t_tr),
        file.path(TAB, "01_ideal.rds"))

cat("\nGATE:", ifelse(all(taus > 0.5), "PASS", "FAIL -- stop and fix before Phase 0 step 3"), "\n")
