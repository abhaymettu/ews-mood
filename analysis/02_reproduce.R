## analysis/02_reproduce.R -- reproduce the published early-warning result
##
## Target: van de Leemput et al. (PNAS 2014) / Wichers et al. (Psychother
## Psychosom 2016) report rising autocorrelation and rising variance in the
## window preceding this participant's transition.
##
## THE INFERENCE PROBLEM. Phase 0 (D9) banned nominal Kendall p-values: on
## overlapping rolling windows a nominally 5% test rejects 28-43% of the time.
## In Phase 0 we replaced them with thresholds calibrated against a simulated
## null. Here there is one person and one transition, so there is no null to
## simulate without assuming the generative model we are trying to test.
##
## Instead: a moving-block bootstrap. Resample blocks of the detrended series
## with replacement, reassign the original observation times, recompute the
## rolling indicator and its tau. This preserves the marginal distribution and
## the within-block dependence -- the local dynamics -- while destroying any
## slow trend in the indicator. It is a null of "these dynamics, in no
## particular order", which is exactly the null the claim needs to beat.
##
## Block length is swept (1 day and 3 days of beeps) because it is a choice.
##
## Data: Kossakowski et al. (2017), JOPD 5:1, doi 10.5334/jopd.29, CC-BY.

source("R/setup.R")
set.seed(SEED)

D <- readRDS(file.path(PROJ, "data", "esm_clean.rds"))
d <- D$data
NBOOT <- 999

## --- moving-block bootstrap null for the trend statistic ---------------------
NCORE <- max(1L, parallel::detectCores() - 1L)
boot_tau <- function(y, t, indicator, win_days, step_days, block_n, nboot = NBOOT) {
  n <- length(y)
  nb <- ceiling(n / block_n)
  starts_all <- seq_len(n - block_n + 1)
  ## 999 surrogates x 4 indicators x 16 series is ~64k rolling passes, and the
  ## OU indicator fits an MLE per window. Parallel by necessity, not taste.
  unlist(parallel::mclapply(seq_len(nboot), function(b) {
    idx <- unlist(lapply(sample(starts_all, nb, replace = TRUE),
                         function(s) s:(s + block_n - 1)))[1:n]
    kendall_tau(rolling(y[idx], t, INDICATORS[[indicator]], win_days, step_days))
  }, mc.cores = NCORE))
}

analyse <- function(y, t, label, win_days = 30, step_days = 2,
                    detrend_method = "gaussian", bw_days = 14, block_n = 6) {
  pre <- detrend(y, t, detrend_method, bw_days)
  out <- lapply(c("variance", "ac1_naive", "ac1_withinday", "ac1_ou"), function(k) {
    r <- rolling(pre, t, INDICATORS[[k]], win_days, step_days)
    if (is.null(r)) return(NULL)
    tau <- kendall_tau(r)
    nul <- boot_tau(pre, t, k, win_days, step_days, block_n)
    data.frame(series = label, indicator = k, tau = tau,
               p_boot = (1 + sum(nul >= tau, na.rm = TRUE)) / (1 + sum(is.finite(nul))),
               null_q95 = quantile(nul, .95, na.rm = TRUE),
               p_nominal = {                       # reported ONLY to show the gap
                 nn <- sum(is.finite(r$val))
                 z <- 3 * tau * sqrt(nn * (nn - 1)) / sqrt(2 * (2 * nn + 5))
                 1 - pnorm(z)
               },
               n_windows = sum(is.finite(r$val)))
  })
  do.call(rbind, out)
}

## --- 1. the headline reproduction -------------------------------------------
## Window: everything up to the transition. Both transition definitions are run
## because D21 found they disagree by 14 days.
res <- list()
for (tdef in names(D$transition)) {
  tday <- D$transition[[tdef]]
  sub <- d[d$day <= tday, ]
  for (v in c("pa", "na_")) {
    r <- analyse(sub[[v]], sub$t_hours, sprintf("%s | pre-transition (%s, day<=%d)",
                                                v, tdef, tday))
    r$transition_def <- tdef; r$series_var <- v; r$n_obs <- nrow(sub)
    res[[length(res) + 1]] <- r
  }
}
res <- do.call(rbind, res)
rownames(res) <- NULL

cat("=== 1. Pre-transition window, composites ===\n")
print(res[c("series_var", "transition_def", "indicator", "n_obs", "n_windows",
            "tau", "p_boot", "p_nominal")], row.names = FALSE, digits = 3)

cat("\nnominal vs bootstrap p-values, same data:\n")
cat(sprintf("  nominal  significant at .05: %d of %d\n", sum(res$p_nominal < .05), nrow(res)))
cat(sprintf("  bootstrap significant at .05: %d of %d\n", sum(res$p_boot < .05), nrow(res)))

## --- 2. per-item: Phase 0's floor/ceiling prediction -------------------------
## PHASE0.md section 4: the variance indicator is confounded by where an item's
## mean sits relative to its response boundaries. Three items here sit within
## 1.5 categories of their floor (suspic 0.26, doubt 0.85, irritat 1.24). If
## Phase 0 is right about the mechanism, those items should behave differently
## from the mid-scale ones. This is a genuine out-of-sample prediction: it was
## made in simulation before this dataset was downloaded.
tday <- D$transition$changepoint
sub <- d[d$day <= tday, ]
items <- do.call(rbind, lapply(D$mood_items, function(m) {
  r <- analyse(sub[[m]], sub$t_hours, m)
  r$item <- m; r
}))
hr <- D$item_ranges
items$headroom <- pmin(hr$headroom_lo, hr$headroom_hi)[match(items$item, hr$item)]
items$near_bound <- items$headroom < 1.5

cat("\n=== 2. Per-item, pre-transition (changepoint definition) ===\n")
iv <- subset(items, indicator == "variance")
print(iv[order(iv$headroom), c("item", "headroom", "tau", "p_boot")],
      row.names = FALSE, digits = 3)

cat("\nmean tau by indicator and scale position:\n")
agg <- aggregate(items["tau"], items[c("indicator", "near_bound")], mean)
print(reshape(agg, idvar = "indicator", timevar = "near_bound", direction = "wide"),
      row.names = FALSE, digits = 3)

## --- 3. full series, for the figure ------------------------------------------
roll_full <- do.call(rbind, lapply(c("pa", "na_"), function(v) {
  pre <- detrend(d[[v]], d$t_hours, "gaussian", 14)
  do.call(rbind, lapply(c("variance", "ac1_naive", "ac1_ou"), function(k) {
    r <- rolling(pre, d$t_hours, INDICATORS[[k]], 30, 2)
    r$indicator <- k; r$series <- v; r
  }))
}))
roll_full$day <- roll_full$end / 24

write.csv(res,   file.path(TAB, "P1_reproduce_composites.csv"), row.names = FALSE)
write.csv(items, file.path(TAB, "P1_reproduce_items.csv"), row.names = FALSE)
saveRDS(list(res = res, items = items, roll = roll_full),
        file.path(TAB, "P1_reproduce.rds"))

## --- figure ------------------------------------------------------------------
lab <- c(variance = "Variance", ac1_naive = "AC1 (naive)", ac1_ou = "AC1 (OU)")
roll_full$indicator <- factor(roll_full$indicator, names(lab), lab)
roll_full$series <- factor(roll_full$series, c("pa", "na_"),
                           c("positive affect", "negative affect"))

p <- ggplot(roll_full, aes(day, val, colour = series)) +
  geom_vline(xintercept = D$transition$changepoint, linetype = 2, colour = "firebrick") +
  geom_vline(xintercept = D$transition$crossing, linetype = 3, colour = "firebrick") +
  geom_rect(data = D$phases[D$phases$phase == 3, ],
            aes(xmin = first, xmax = last, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey70", alpha = .25) +
  geom_line(linewidth = .5) +
  facet_wrap(~indicator, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("positive affect" = "steelblue4",
                                 "negative affect" = "darkorange3")) +
  labs(title = "Early warning indicators, Kossakowski et al. (2017) participant",
       subtitle = "30-day rolling windows. Grey band = double-blind dose reduction. Dashed = changepoint (day 127), dotted = mean+2SD crossing (day 141).",
       x = "day of study", y = NULL, colour = NULL,
       caption = "Data: Kossakowski et al. 2017, JOPD 5:1, CC-BY 4.0") +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "P1_indicators.png"), p, width = 10, height = 8, dpi = 150)

sclp <- ggplot(D$scl, aes(day, dep)) +
  geom_vline(xintercept = D$transition$changepoint, linetype = 2, colour = "firebrick") +
  geom_vline(xintercept = D$transition$crossing, linetype = 3, colour = "firebrick") +
  geom_line(colour = "grey40") + geom_point(size = 1.2) +
  labs(title = "Weekly SCL-90 depression score",
       subtitle = "The rise is gradual. Calling it a critical transition is an interpretation (D21).",
       x = "day of study", y = "SCL-90 depression")
ggsave(file.path(FIG, "P1_scl90.png"), sclp, width = 9, height = 3.5, dpi = 150)
