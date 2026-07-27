## analysis/06_replicate_artifact.R -- does the discretisation artefact replicate?
##
## Phase 1 found, in one participant, that the rolling-variance early warning
## signal was concentrated in items sitting near their response floor, and was
## not distinguishable from what rounding alone produces. One participant is not
## evidence that this happens generally.
##
## Fisher et al. (2017) provides the ideal test, for a reason that has nothing to
## do with transitions: **their items are 0-100 visual analog scales.** With 101
## categories the rounding artefact should be negligible. Discretising those same
## series to 7 categories should create it.
##
## That is a within-dataset, within-participant, within-item comparison: identical
## underlying measurements, two response resolutions, nothing else changed. No
## transition is needed, because the claim being tested is about the measurement
## model, not about critical transitions.
##
## PREDICTION (from PHASE0.md section 4, made before this dataset was obtained):
##   on 0-100 VAS      -- variance tau should be unrelated to an item's position
##                        on the scale
##   on 7-point Likert -- variance tau should be higher for items sitting near a
##                        scale bound, reproducing the Kossakowski pattern
##
## Data: Fisher, A. J., Reeves, J. W., Lawyer, G., Medaglia, J. D., & Rubel, J. A.
## (2017). Exploring the idiographic dynamics of mood and anxiety via network
## analysis. Journal of Abnormal Psychology, 126(8), 1044-1056. OSF: osf.io/zefbc

source("R/setup.R")
set.seed(SEED)

FILES <- list.files(file.path(PROJ, "data-raw", "fisher", "R Data"),
                    full.names = TRUE, pattern = "[.]RData$")
cat("participants:", length(FILES), "\n")

WIN_H <- 10 * 24      # 10-day rolling window
STEP_H <- 24          # 1-day step
BEEP_H <- 6           # ~4 prompts/day

to_likert <- function(y, K = 7) round((y - 0) / 100 * (K - 1)) + 1

one_participant <- function(f) {
  e <- new.env(); load(f, envir = e)
  if (!exists("trim", envir = e)) return(NULL)
  d <- get("trim", envir = e)
  if (nrow(d) < 80) return(NULL)
  t <- (seq_len(nrow(d)) - 1) * BEEP_H
  pid <- sub("_final[.]RData$", "", basename(f))

  do.call(rbind, lapply(names(d), function(item) {
    y <- d[[item]]
    if (!is.numeric(y)) return(NULL)
    ## listwise on this item only: the detrend smoother cannot carry NAs, and
    ## dropping them keeps the true clock via the retained timestamps
    ok <- is.finite(y)
    if (sum(ok) < 80 || sd(y[ok]) == 0) return(NULL)
    yy <- y[ok]; tt <- t[ok]
    yl <- to_likert(yy)

    tau_of <- function(v) kendall_tau(rolling(detrend(v, tt, "gaussian", 5), tt,
                                              ind_var, WIN_H / 24, STEP_H / 24))
    tv <- tau_of(yy); tl <- tau_of(yl)
    if (!is.finite(tv) || !is.finite(tl)) return(NULL)

    ## position of the item's mean on the 7-point scale
    m <- mean(yl)
    data.frame(pid = pid, item = item,
               mean_vas = mean(yy),
               headroom = min(m - 1, 7 - m),          # categories to nearest bound
               tau_vas = tv, tau_likert = tl, d_tau = tl - tv)
  }))
}

res <- do.call(rbind, lapply(FILES, one_participant))
cat(sprintf("usable item-series: %d across %d participants\n",
            nrow(res), length(unique(res$pid))))

res$near_bound <- res$headroom < 1.5

## --- the registered comparison ----------------------------------------------
cat("\n=== Mean variance-trend tau by scale position and response resolution ===\n")
agg <- aggregate(res[c("tau_vas", "tau_likert")], list(near_bound = res$near_bound), mean)
agg$n <- as.vector(table(res$near_bound))
names(agg)[1] <- "within 1.5 categories of a bound"
print(agg, row.names = FALSE, digits = 3)

sp_vas <- diff(rev(agg$tau_vas)); sp_lik <- diff(rev(agg$tau_likert))
cat(sprintf("\nsplit (near-bound minus mid-scale):  VAS %+.3f   Likert-7 %+.3f\n",
            -sp_vas, -sp_lik))

## Does discretisation systematically inflate tau for boundary-adjacent items?
cat("\n=== Effect of discretising the SAME series ===\n")
m <- lm(d_tau ~ headroom, data = res)
print(summary(m)$coefficients, digits = 3)
ct_v <- cor.test(res$tau_vas, res$headroom, method = "kendall")
ct_l <- cor.test(res$tau_likert, res$headroom, method = "kendall")
cat(sprintf("\ncorrelation of variance-tau with headroom:\n  VAS      tau=%+.3f p=%.4f\n  Likert-7 tau=%+.3f p=%.4f\n",
            ct_v$estimate, ct_v$p.value, ct_l$estimate, ct_l$p.value))

## per-participant, so the result is not driven by a few people
pp <- do.call(rbind, lapply(split(res, res$pid), function(z) {
  if (length(unique(z$near_bound)) < 2) return(NULL)
  data.frame(pid = z$pid[1],
             vas = mean(z$tau_vas[z$near_bound]) - mean(z$tau_vas[!z$near_bound]),
             lik = mean(z$tau_likert[z$near_bound]) - mean(z$tau_likert[!z$near_bound]))
}))
cat(sprintf("\nparticipants where the near-bound split is larger under Likert: %d of %d\n",
            sum(pp$lik > pp$vas), nrow(pp)))
cat(sprintf("paired test on the per-participant split (Likert vs VAS): p = %.5f\n",
            wilcox.test(pp$lik, pp$vas, paired = TRUE)$p.value))

write.csv(res, file.path(TAB, "P4_fisher_replication.csv"), row.names = FALSE)
saveRDS(list(res = res, pp = pp, agg = agg), file.path(TAB, "P4_fisher.rds"))

## --- figure ------------------------------------------------------------------
long <- rbind(
  data.frame(headroom = res$headroom, tau = res$tau_vas, scale = "0-100 VAS (as collected)"),
  data.frame(headroom = res$headroom, tau = res$tau_likert, scale = "same data, 7-point Likert"))
p <- ggplot(long, aes(headroom, tau)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(alpha = .25, size = .8) +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = .8) +
  facet_wrap(~scale) +
  labs(title = "Discretising the same measurements creates the variance artefact",
       subtitle = sprintf("%d item-series from %d participants (Fisher et al. 2017). Left: as collected. Right: identical data rounded to 7 categories.",
                          nrow(res), length(unique(res$pid))),
       x = "categories from the item mean to the nearest scale bound",
       y = "Kendall tau of rolling variance",
       caption = "Data: Fisher et al. 2017, osf.io/zefbc")
ggsave(file.path(FIG, "P4_fisher_replication.png"), p, width = 10, height = 4.8, dpi = 150)
