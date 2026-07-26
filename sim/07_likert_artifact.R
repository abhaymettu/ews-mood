## sim/07_likert_artifact.R -- the response scale can invent, or reverse, the
## variance early-warning signal.
##
## Found while checking why the clinical (drift) null had a NEGATIVE variance
## trend (tau = -0.29) when the latent variance in those same runs rises by 14%,
## exactly as critical slowing down predicts.
##
## Mechanism. The variance of a rounded variable depends on where its mean sits
## between two category boundaries. A mean sitting exactly on a boundary
## produces near-Bernoulli variance across the two categories; a mean sitting in
## the middle of a category produces almost none. Approaching a transition the
## latent mean *slides*, so it moves relative to the boundaries, and the rounding
## variance moves with it -- by an amount that can dwarf the real change in
## process variance.
##
## On a single latent path, moving the response scale alone drove the variance
## Kendall tau from -0.98 to +0.96.
##
## This matters because variance was the strongest indicator in the power grid
## (sim/04). If part of that strength is this artefact, the grid's headline
## recommendation is contaminated. That is what this script tests: does the
## variance advantage survive when discretisation is removed?
##
## Metric is AUC of transition vs the clinical null, which is threshold-free.

source("R/setup.R")

NREP <- 200
DAYS <- 120; PPD <- 6; SM <- 0.1
INDS <- c("variance", "ac1_naive", "ac1_withinday", "ac1_ou")

## K = 0 is continuous (no rounding). centre 4 / scale 2 puts the euthymic state
## at ~6.4 on a 1-7 item, i.e. hard against the 6/7 boundary and near ceiling --
## the placement used in sim/04, chosen arbitrarily. centre 3.5 / scale 1.3 puts
## euthymic ~5.0 and depressed ~2.0, using the middle of the scale.
CFG <- list(
  list(id = "continuous",           K = 0,  scale_pts = 2.0, center = 4.0),
  list(id = "Likert-7, at ceiling", K = 7,  scale_pts = 2.0, center = 4.0),
  list(id = "Likert-7, centred",    K = 7,  scale_pts = 1.3, center = 3.5),
  list(id = "Likert-11, centred",   K = 11, scale_pts = 2.2, center = 6.0)
)

run <- function(i, scenario) {
  set.seed(SEED + 100000 * match(scenario, c("transition", "static", "drift")) + i)
  p <- simulate_latent(SIM_DAYS, scenario); tt <- find_transition(p)
  if (scenario == "transition") {
    if (is.na(tt)) return(NULL); end_h <- tt
  } else {
    if (!is.na(tt) && tt <= REF_DAY * 24) return(NULL); end_h <- REF_DAY * 24
  }
  do.call(rbind, lapply(CFG, function(cf) {
    set.seed(SEED + 900000 + i)
    o <- observe(p, end_h, DAYS, PPD, sigma_m = SM, K = cf$K, scale_pts = cf$scale_pts)
    if (is.null(o)) return(NULL)
    ## observe() uses centre 4; re-render the response with this config's centre
    set.seed(SEED + 900000 + i)
    o$y <- respond(o$x, sigma_m = SM, scale_pts = cf$scale_pts,
                   center = cf$center, K = cf$K)
    pre <- detrend(o$y, o$t, "gaussian", bw_days = DAYS / 8)
    data.frame(rep = i, scenario = scenario, cfg = cf$id, indicator = INDS,
               tau = vapply(INDS, function(k)
                 ews_tau(NULL, o$t, k, win_days = DAYS / 2, step_days = 4, pre = pre),
                 numeric(1)),
               ceiling_frac = mean(o$y >= max(cf$K, 1) & cf$K > 0))
  }))
}

nc <- max(1L, parallel::detectCores() - 1L)
res <- do.call(rbind, unlist(lapply(c("transition", "static", "drift"), function(s)
  parallel::mclapply(1:NREP, run, scenario = s, mc.cores = nc)), recursive = FALSE))

perf <- do.call(rbind, lapply(split(res, res[c("cfg", "indicator")], drop = TRUE), function(z) {
  tr <- na.omit(z$tau[z$scenario == "transition"])
  st <- na.omit(z$tau[z$scenario == "static"])
  dr <- na.omit(z$tau[z$scenario == "drift"])
  if (length(tr) < 30 || length(st) < 30) return(NULL)
  thr <- quantile(st, 0.95)
  data.frame(cfg = z$cfg[1], indicator = z$indicator[1],
             tau_transition = mean(tr), tau_drift = mean(dr),
             power = mean(tr > thr), fpr_drift = mean(dr > thr),
             auc = mean(outer(tr, dr, ">")) + 0.5 * mean(outer(tr, dr, "==")))
}))
rownames(perf) <- NULL
perf <- perf[order(perf$cfg, -perf$auc), ]
cat("=== Does the variance advantage survive removing discretisation? ===\n")
print(perf, row.names = FALSE, digits = 3)

cat("\n=== mean tau on the CLINICAL NULL (should be mildly positive: it really does slow) ===\n")
dn <- perf[perf$indicator == "variance", c("cfg", "tau_drift")]
print(dn, row.names = FALSE, digits = 3)
cat("\nA negative value here is the artefact: the latent variance in these runs rises ~14%.\n")

write.csv(perf, file.path(TAB, "07_likert_artifact.csv"), row.names = FALSE)
saveRDS(res, file.path(TAB, "07_raw.rds"))

p <- ggplot(perf, aes(cfg, auc, fill = indicator)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0.3, 1)) +
  labs(title = "The response scale changes which indicator looks best",
       subtitle = "AUC, transition vs a deteriorating-but-not-tipping null; 0.5 = useless",
       x = NULL, y = "AUC", fill = NULL) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(FIG, "07_likert_artifact.png"), p, width = 9, height = 5, dpi = 150)
