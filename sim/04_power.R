## sim/04_power.R -- Phase 0 step 4: detection power and false-alarm rate as a
## function of the things an ESM designer actually controls.
##
## DESIGN
## ------
## Three scenarios, all run through the identical observation grid:
##   transition -- tips. Observation window = the D days ENDING at the tipping
##                 point. That is the situation a clinician is in: D days of
##                 history and an event at the end.
##   static     -- c constant. Used ONLY to set the decision threshold, per cell.
##   drift      -- c descends but never reaches the fold. The realistic false
##                 alarm: a person who deteriorates and does not become ill.
##
## Thresholds are set per cell from the static null, so every cell is a genuine
## 5%-level test and the power numbers are comparable across cells. Using one
## global threshold, or a nominal p-value, would make short/sparse cells look
## powerful purely because their tau is noisier (sim/02).
##
## The latent path is simulated ONCE per replicate and every observation
## configuration is derived from it. Conditions are therefore paired on the
## underlying trajectory, which removes the between-path variance (SD ~ 0.33,
## see sim/03) from all within-replicate comparisons.

source("R/setup.R")

NREP     <- 400

GRID <- expand.grid(
  days      = c(30, 60, 120, 180),
  ppd       = c(3, 6, 10),
  miss      = c("none", "mcar30", "mnar30"),
  sigma_m   = c(0.1, 0.4, 0.8),
  stringsAsFactors = FALSE
)
WIN_FRAC <- c(0.25, 0.50)   # rolling window as a fraction of series length
INDS     <- c("variance", "ac1_naive", "ac1_withinday", "ac1_ou")
MISS     <- list(none = c("none", 0), mcar30 = c("mcar", 0.30), mnar30 = c("mnar", 0.30))

## One replicate: simulate a path, evaluate every grid cell on it.
run_rep <- function(i, scenario) {
  set.seed(SEED + 100000 * match(scenario, c("transition", "static", "drift")) + i)
  p <- simulate_latent(SIM_DAYS, scenario)
  tt <- find_transition(p)
  if (scenario == "transition") {
    if (is.na(tt)) return(NULL)
    end_h <- tt
  } else {
    if (!is.na(tt) && tt <= REF_DAY * 24) return(NULL)   # null that tipped: discard
    end_h <- REF_DAY * 24
  }

  ## preallocated columns: rbind-ing ~900 one-row frames per replicate costs
  ## more than the statistics do
  N <- nrow(GRID) * length(WIN_FRAC) * length(INDS)
  gi <- integer(N); wf_ <- numeric(N); ind_ <- character(N)
  tau_ <- rep(NA_real_, N); nob_ <- integer(N); j <- 0L

  for (g in seq_len(nrow(GRID))) {
    G <- GRID[g, ]
    if (end_h < G$days * 24) next                        # not enough history
    m <- MISS[[G$miss]]
    set.seed(SEED + 900000 + i * 1000 + g)               # observation-layer noise
    o <- observe(p, end_h, G$days, G$ppd, sigma_m = G$sigma_m, K = 7,
                 mechanism = m[1], miss_rate = as.numeric(m[2]))
    if (is.null(o) || nrow(o) < 60) next
    ## detrend bandwidth scales with series length: a fixed 14-day bandwidth
    ## removes almost everything from a 30-day series and almost nothing from a
    ## 180-day one, which would confound length with detrending strength (D4).
    pre <- detrend(o$y, o$t, "gaussian", bw_days = max(7, G$days / 8))
    for (wf in WIN_FRAC) {
      w <- G$days * wf
      for (k in INDS) {
        j <- j + 1L
        gi[j] <- g; wf_[j] <- wf; ind_[j] <- k; nob_[j] <- nrow(o)
        tau_[j] <- ews_tau(NULL, o$t, k, win_days = w,
                           step_days = max(1, w / 15), pre = pre)
      }
    }
  }
  if (!j) return(NULL)
  s <- seq_len(j)
  cbind(GRID[gi[s], ], win_frac = wf_[s], win_days = GRID$days[gi[s]] * wf_[s],
        indicator = ind_[s], rep = i, scenario = scenario, n_obs = nob_[s],
        tau = tau_[s])
}

ncore <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("grid: %d cells x %d window widths x %d indicators; %d reps x 3 scenarios on %d cores\n",
            nrow(GRID), length(WIN_FRAC), length(INDS), NREP, ncore))

t0 <- Sys.time()
all_res <- do.call(rbind, unlist(lapply(c("transition", "static", "drift"), function(s) {
  cat("  scenario:", s, "\n")
  parallel::mclapply(1:NREP, run_rep, scenario = s, mc.cores = ncore)
}), recursive = FALSE))
cat(sprintf("done in %.1f min; %d rows\n",
            as.numeric(Sys.time() - t0, units = "mins"), nrow(all_res)))
saveRDS(all_res, file.path(TAB, "04_power_raw.rds"))

## --- per-cell threshold from the static null, then power and false alarms ----
key <- c("days", "ppd", "miss", "sigma_m", "win_frac", "indicator")
sp  <- split(all_res, all_res[key], drop = TRUE)

power <- do.call(rbind, lapply(sp, function(z) {
  tr <- z$tau[z$scenario == "transition"]
  st <- z$tau[z$scenario == "static"]
  dr <- z$tau[z$scenario == "drift"]
  if (sum(is.finite(st)) < 50 || sum(is.finite(tr)) < 50) return(NULL)
  thr <- quantile(st, 0.95, na.rm = TRUE)
  data.frame(z[1, key],
             n_obs   = round(mean(z$n_obs)),
             thr     = as.numeric(thr),
             power   = mean(tr > thr, na.rm = TRUE),
             fpr_drift = mean(dr > thr, na.rm = TRUE),
             ## AUC vs the clinical null: threshold-free, so it is not hostage to
             ## the 5% convention. 0.5 = useless.
             auc     = {
               a <- na.omit(tr); b <- na.omit(dr)
               if (!length(a) || !length(b)) NA_real_ else
                 mean(outer(a, b, ">")) + 0.5 * mean(outer(a, b, "=="))
             },
             n_tr = sum(is.finite(tr)), n_null = sum(is.finite(st)))
}))
rownames(power) <- NULL
write.csv(power, file.path(TAB, "04_power.csv"), row.names = FALSE)

cat("\n=== POWER at a per-cell 5% level, best indicator per configuration ===\n")
best <- do.call(rbind, lapply(split(power, power[c("days", "ppd", "miss", "sigma_m")], drop = TRUE),
                              function(z) z[which.max(z$auc), ]))
print(best[order(-best$power),
           c("days", "ppd", "miss", "sigma_m", "indicator", "win_frac",
             "n_obs", "power", "fpr_drift", "auc")][1:20, ], row.names = FALSE, digits = 3)

cat("\n=== marginal effects (mean power over the rest of the grid) ===\n")
for (v in c("days", "ppd", "miss", "sigma_m", "indicator", "win_frac")) {
  a <- aggregate(cbind(power, fpr_drift, auc) ~ get(v), power, mean)
  names(a)[1] <- v
  cat("\n"); print(a, row.names = FALSE, digits = 3)
}

## --- figures -----------------------------------------------------------------
power$miss <- factor(power$miss, c("none", "mcar30", "mnar30"))
pw <- subset(power, win_frac == 0.5)

p1 <- ggplot(pw, aes(days, power, colour = factor(ppd), group = ppd)) +
  geom_hline(yintercept = 0.05, linetype = 3, colour = "grey50") +
  geom_hline(yintercept = 0.80, linetype = 2, colour = "grey50") +
  geom_line() + geom_point(size = 1) +
  facet_grid(sigma_m ~ miss, labeller = label_both) +
  scale_colour_viridis_d(end = .85, name = "prompts/day") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Power to detect critical slowing down before a transition",
       subtitle = "5% level calibrated per cell against a no-change null; best-case indicator shown separately below",
       x = "days of history before the event", y = "power") +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "04_power_grid.png"), p1, width = 10, height = 7.5, dpi = 150)

p2 <- ggplot(pw, aes(days, power, colour = indicator)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey50") +
  geom_line() +
  facet_grid(sigma_m ~ ppd, labeller = label_both) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Indicator comparison", x = "days of history", y = "power") +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "04_power_indicators.png"), p2, width = 10, height = 7, dpi = 150)

p3 <- ggplot(pw, aes(fpr_drift, power, colour = factor(days), shape = miss)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70") +
  geom_point(alpha = .8) +
  scale_colour_viridis_d(end = .85, name = "days") +
  labs(title = "Power against realistic false alarms",
       subtitle = "x-axis: alarm rate on a person who deteriorates but never tips. Points near the diagonal are useless.",
       x = "false alarm rate (clinical null)", y = "power")
ggsave(file.path(FIG, "04_power_vs_falsealarm.png"), p3, width = 8, height = 5.5, dpi = 150)
