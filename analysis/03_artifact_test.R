## analysis/03_artifact_test.R -- is the rising variance critical slowing down,
## or is it the mean sliding across response-category boundaries?
##
## THE PROBLEM. analysis/02 found that the variance early-warning signal in this
## participant is concentrated entirely in three items -- suspic, doubt,
## irritat -- which are exactly the three sitting within 1.5 categories of their
## floor. Mean tau for the variance indicator was +0.62 for floored items and
## -0.07 for mid-scale items. The autocorrelation indicators show no such split.
##
## That matches the artefact PHASE0.md section 4 predicted in simulation, before
## this dataset was downloaded. But it is ALSO exactly what genuine symptom
## emergence would look like: those three items are floored because the person
## is not yet symptomatic, and as they become depressed the items lift off the
## floor and genuinely become more variable. Both stories predict the same
## correlation. The correlation alone proves nothing.
##
## THE TEST. Take each item's observed rolling mean trajectory. Generate a
## synthetic series that follows that same mean trajectory but has CONSTANT
## latent variance by construction -- no critical slowing down, none. Round it
## to the item's integer response scale. Then run the identical detrend-and-roll
## pipeline and see how much variance trend appears.
##
## Any tau in the synthetic series is pure rounding artefact, because there is
## no real variance change in it to detect.
##
##   synthetic tau ~ observed tau  -> the artefact is sufficient; no slowing
##                                    down needs to be invoked
##   synthetic tau << observed tau -> something real is happening on top of it
##
## Data: Kossakowski et al. (2017), JOPD 5:1, doi 10.5334/jopd.29, CC-BY.

source("R/setup.R")
set.seed(SEED)

D <- readRDS(file.path(PROJ, "data", "esm_clean.rds"))
d <- D$data
NSIM <- 300
WIN <- 30; STEP <- 2; BW <- 14

tday <- D$transition$changepoint
sub  <- d[d$day <= tday, ]
t    <- sub$t_hours

## Gaussian smoother fit -- the thing detrend() subtracts. This IS the item's
## mean trajectory, and it is what we hold onto while discarding everything else.
smooth_fit <- function(y, t, bw_days = BW) y - detrend(y, t, "gaussian", bw_days)

## Item response bounds, from the codebook anchoring recorded in 01_clean.
bounds <- function(item) {
  lo <- if (item %in% c("mood_down", "mood_lonely", "mood_anxious", "mood_guilty")) -3 else 1
  c(lo, lo + 6)
}

res <- do.call(rbind, lapply(D$mood_items, function(item) {
  y  <- sub[[item]]
  mu <- smooth_fit(y, t)                       # observed mean trajectory
  s  <- sd(y - mu)                             # residual SD, held CONSTANT below
  b  <- bounds(item)

  obs_tau <- kendall_tau(rolling(detrend(y, t, "gaussian", BW), t, ind_var, WIN, STEP))

  ## synthetic: same mean path, constant latent variance, same rounding
  sim_tau <- unlist(parallel::mclapply(seq_len(NSIM), function(i) {
    ys <- pmin(pmax(round(mu + rnorm(length(mu), 0, s)), b[1]), b[2])
    kendall_tau(rolling(detrend(ys, t, "gaussian", BW), t, ind_var, WIN, STEP))
  }, mc.cores = max(1L, parallel::detectCores() - 1L)))

  hr <- D$item_ranges
  data.frame(item = item,
             headroom = pmin(hr$headroom_lo, hr$headroom_hi)[match(item, hr$item)],
             mean_shift = max(mu) - min(mu),
             obs_tau = obs_tau,
             sim_tau_mean = mean(sim_tau, na.rm = TRUE),
             sim_tau_q05 = quantile(sim_tau, .05, na.rm = TRUE),
             sim_tau_q95 = quantile(sim_tau, .95, na.rm = TRUE),
             ## does the constant-variance synthetic already reach the observed?
             p_artefact_enough = mean(sim_tau >= obs_tau, na.rm = TRUE))
}))
rownames(res) <- NULL
res <- res[order(res$headroom), ]

cat("=== Can a CONSTANT-variance series reproduce the observed variance trend? ===\n")
cat("(sim_tau is generated with no change in latent variance whatsoever)\n\n")
print(res[c("item", "headroom", "mean_shift", "obs_tau", "sim_tau_mean",
            "sim_tau_q95", "p_artefact_enough")], row.names = FALSE, digits = 3)

near <- res$headroom < 1.5
cat(sprintf("\nfloored items (n=%d):   observed tau %+.3f   artefact-only tau %+.3f  -> %.0f%% explained\n",
            sum(near), mean(res$obs_tau[near]), mean(res$sim_tau_mean[near]),
            100 * mean(res$sim_tau_mean[near]) / mean(res$obs_tau[near])))
cat(sprintf("mid-scale items (n=%d): observed tau %+.3f   artefact-only tau %+.3f\n",
            sum(!near), mean(res$obs_tau[!near]), mean(res$sim_tau_mean[!near])))
cat(sprintf("\nitems whose observed trend is NOT beyond the artefact (p > .05): %d of %d\n",
            sum(res$p_artefact_enough > .05), nrow(res)))

write.csv(res, file.path(TAB, "P1_artefact_test.csv"), row.names = FALSE)
saveRDS(res, file.path(TAB, "P1_artefact_test.rds"))

res$grp <- ifelse(near, "within 1.5 categories of floor", "mid-scale")
p <- ggplot(res, aes(reorder(item, headroom))) +
  geom_linerange(aes(ymin = sim_tau_q05, ymax = sim_tau_q95), colour = "grey60", linewidth = 3) +
  geom_point(aes(y = sim_tau_mean), colour = "grey30", size = 2) +
  geom_point(aes(y = obs_tau, colour = grp), size = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  coord_flip() +
  scale_colour_manual(values = c("within 1.5 categories of floor" = "firebrick",
                                 "mid-scale" = "steelblue4")) +
  labs(title = "Observed variance trend vs what rounding alone produces",
       subtitle = "Grey = 90% range of tau from a synthetic series with CONSTANT latent variance.\nWhere the coloured point sits inside the grey band, no slowing down needs to be invoked.",
       x = NULL, y = "Kendall tau of rolling variance", colour = NULL,
       caption = "Data: Kossakowski et al. 2017, JOPD 5:1, CC-BY 4.0") +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "P1_artefact_test.png"), p, width = 9, height = 5.5, dpi = 150)
