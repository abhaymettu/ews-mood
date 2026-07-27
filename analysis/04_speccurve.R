## analysis/04_speccurve.R -- Phase 2 specification curve
##
## Preregistered in preregistration/PREREG_phase2_phase3.md, committed before
## this file existed. Every dimension and level below is as registered; any
## departure goes in preregistration/DEVIATIONS.md.
##
## 4 windows x 5 detrends x 4 series x 3 spacings x 3 taper treatments
## x 2 transition definitions x 4 indicators = 5,760 specifications.
##
## Inference is the registered joint test: 100 moving-block-bootstrap surrogate
## datasets, the ENTIRE curve recomputed on each, giving a null distribution for
## the curve's median tau (S1) and proportion positive (S2). Nominal p-values
## appear nowhere (D9).
##
## Data: Kossakowski et al. (2017), JOPD 5:1, doi 10.5334/jopd.29, CC-BY.

source("R/setup.R")
set.seed(SEED)

D <- readRDS(file.path(PROJ, "data", "esm_clean.rds"))
d <- D$data

FLOORED  <- c("mood_suspic", "mood_doubt", "mood_irritat")
MIDSCALE <- setdiff(D$mood_items, FLOORED)

SPEC <- expand.grid(
  win        = c(20, 30, 45, 60),
  detrend    = c("none", "linear", "gauss7", "gauss14", "gauss28"),
  series     = c("pa", "na_", "floored", "midscale"),
  spacing    = c("native", "withinday", "interp3h"),
  taper      = c("ignored", "covariate", "break"),
  transition = c("changepoint", "crossing"),
  indicator  = c("variance", "ac1_naive", "ac1_withinday", "ac1_ou"),
  stringsAsFactors = FALSE
)
cat(sprintf("specifications: %d\n", nrow(SPEC)))

## --- build each series once --------------------------------------------------
mk_series <- function(dat, which) switch(which,
  pa       = dat$pa,
  na_      = dat$na_,
  floored  = rowMeans(dat[FLOORED]),
  midscale = rowMeans(dat[MIDSCALE])
)

## taper treatments (D20). "covariate" residualises on the taper indicators;
## "break" centres each study phase separately, which is what treating the taper
## as a structural break amounts to for a variance/AC analysis.
apply_taper <- function(y, dat, how) switch(how,
  ignored   = y,
  covariate = residuals(lm(y ~ dat$tapering + dat$post_taper)),
  `break`   = y - ave(y, dat$phase)
)

## spacing treatments (D18). "interp3h" puts the series on a regular 3-hour grid
## by linear interpolation, which is what a lot of applied work silently does.
apply_spacing <- function(y, t, how) switch(how,
  native    = list(y = y, t = t),
  withinday = list(y = y, t = t),         # handled by the indicator itself
  interp3h  = { g <- seq(min(t), max(t), by = 3)
                list(y = approx(t, y, xout = g, rule = 2)$y, t = g) }
)

DTR <- list(none = c("none", 0), linear = c("linear", 0),
            gauss7 = c("gaussian", 7), gauss14 = c("gaussian", 14),
            gauss28 = c("gaussian", 28))

## One specification -> one tau. `dat` is passed in so surrogates can reuse it.
run_spec <- function(k, dat) {
  S <- SPEC[k, ]
  tday <- D$transition[[S$transition]]
  sub  <- dat[dat$day <= tday, ]
  if (nrow(sub) < 100) return(NA_real_)

  y <- mk_series(sub, S$series)
  y <- apply_taper(y, sub, S$taper)
  sp <- apply_spacing(y, sub$t_hours, S$spacing)

  dd <- DTR[[S$detrend]]
  pre <- detrend(sp$y, sp$t, dd[1], as.numeric(dd[2]))

  ind <- if (S$spacing == "withinday" && S$indicator == "ac1_naive")
           "ac1_withinday" else S$indicator
  kendall_tau(rolling(pre, sp$t, INDICATORS[[ind]], S$win, max(1, S$win / 15)))
}

NCORE <- max(1L, parallel::detectCores() - 1L)
curve_of <- function(dat) unlist(parallel::mclapply(seq_len(nrow(SPEC)),
                                                    run_spec, dat = dat,
                                                    mc.cores = NCORE))

## --- observed curve ----------------------------------------------------------
t0 <- Sys.time()
SPEC$tau <- curve_of(d)
cat(sprintf("observed curve in %.1f min; %d of %d specs failed to compute\n",
            as.numeric(Sys.time() - t0, units = "mins"),
            sum(is.na(SPEC$tau)), nrow(SPEC)))

S1_obs <- median(SPEC$tau, na.rm = TRUE)
S2_obs <- mean(SPEC$tau > 0, na.rm = TRUE)
cat(sprintf("\nS1 median tau = %+.3f\nS2 proportion positive = %.3f\n", S1_obs, S2_obs))

## --- registered joint test ---------------------------------------------------
## Resample blocks of the observed series, reassign the original times, and
## recompute the whole curve. 100 surrogates, block = one day of beeps.
NPERM <- 100
block_surrogate <- function(dat, block_n) {
  n <- nrow(dat); nb <- ceiling(n / block_n)
  starts <- sample(seq_len(n - block_n + 1), nb, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) s:(s + block_n - 1)))[1:n]
  out <- dat
  vars <- c(D$mood_items, "pa", "na_")
  out[vars] <- dat[idx, vars]                 # shuffle values, keep the clock
  out
}

joint <- function(block_n, nperm = NPERM) {
  t0 <- Sys.time()
  m <- t(vapply(seq_len(nperm), function(b) {
    set.seed(SEED + 7000 + b)
    tt <- curve_of(block_surrogate(d, block_n))
    c(S1 = median(tt, na.rm = TRUE), S2 = mean(tt > 0, na.rm = TRUE))
  }, numeric(2)))
  cat(sprintf("  block=%d done in %.1f min\n", block_n,
              as.numeric(Sys.time() - t0, units = "mins")))
  data.frame(block_n = block_n,
             S1_null_med = median(m[, "S1"]), S1_null_q95 = quantile(m[, "S1"], .95),
             S2_null_med = median(m[, "S2"]), S2_null_q95 = quantile(m[, "S2"], .95),
             p_S1 = (1 + sum(m[, "S1"] >= S1_obs)) / (1 + nperm),
             p_S2 = (1 + sum(m[, "S2"] >= S2_obs)) / (1 + nperm))
}

cat("\nrunning joint test (100 surrogate curves per block length)...\n")
jt <- rbind(joint(6), joint(18))
cat("\n=== REGISTERED JOINT TEST ===\n")
cat(sprintf("observed  S1 = %+.3f   S2 = %.3f\n", S1_obs, S2_obs))
print(jt, row.names = FALSE, digits = 3)

## --- registered predictions --------------------------------------------------
cat("\n=== REGISTERED PREDICTIONS ===\n")
p1 <- S1_obs > 0 && min(jt$p_S1) > .05
cat(sprintf("P1 median tau positive but p > .05:              %s (tau %+.3f, p %.3f)\n",
            ifelse(p1, "SUPPORTED", "FALSIFIED"), S1_obs, min(jt$p_S1)))

top <- SPEC[!is.na(SPEC$tau), ]
top <- top[top$tau >= quantile(top$tau, .90), ]
frac <- mean(top$indicator == "variance" & top$series == "floored")
p2 <- frac >= 0.5
cat(sprintf("P2 top decile is variance-on-floored-items:      %s (%.0f%%)\n",
            ifelse(p2, "SUPPORTED", "FALSIFIED"), 100 * frac))

rng_of <- function(v) { m <- tapply(SPEC$tau, SPEC[[v]], median, na.rm = TRUE)
                        diff(range(m, na.rm = TRUE)) }
r_dt <- rng_of("detrend"); r_wn <- rng_of("win")
p3 <- r_dt > r_wn
cat(sprintf("P3 detrend moves tau more than window:          %s (%.3f vs %.3f)\n",
            ifelse(p3, "SUPPORTED", "FALSIFIED"), r_dt, r_wn))

flips <- vapply(c("win", "detrend", "series", "spacing", "taper", "transition", "indicator"),
                function(v) { m <- tapply(SPEC$tau, SPEC[[v]], median, na.rm = TRUE)
                              any(m > 0, na.rm = TRUE) && any(m < 0, na.rm = TRUE) }, logical(1))
p4 <- !any(flips[setdiff(names(flips), "series")])
cat(sprintf("P4 only `series` flips the median sign:         %s (flippers: %s)\n",
            ifelse(p4, "SUPPORTED", "FALSIFIED"),
            paste(names(flips)[flips], collapse = ", ")))

## --- what each dimension does ------------------------------------------------
cat("\n=== median tau by level of each dimension ===\n")
for (v in c("series", "indicator", "detrend", "win", "spacing", "taper", "transition")) {
  m <- tapply(SPEC$tau, SPEC[[v]], median, na.rm = TRUE)
  cat(sprintf("%-11s %s\n", v,
              paste(sprintf("%s=%+.3f", names(m), m), collapse = "  ")))
}

write.csv(SPEC, file.path(TAB, "P2_speccurve.csv"), row.names = FALSE)
saveRDS(list(spec = SPEC, joint = jt, S1 = S1_obs, S2 = S2_obs),
        file.path(TAB, "P2_speccurve.rds"))

## --- figure ------------------------------------------------------------------
sc <- SPEC[!is.na(SPEC$tau), ]
sc <- sc[order(sc$tau), ]; sc$rank <- seq_len(nrow(sc))
sc$hl <- ifelse(sc$indicator == "variance" & sc$series == "floored",
                "variance on floored items", "everything else")

p_top <- ggplot(sc, aes(rank, tau, colour = hl)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_hline(yintercept = S1_obs, linetype = 2, colour = "firebrick") +
  geom_point(size = .35) +
  scale_colour_manual(values = c("variance on floored items" = "firebrick",
                                 "everything else" = "grey55")) +
  labs(title = sprintf("Specification curve: %d analytic choices, one participant", nrow(sc)),
       subtitle = sprintf("median tau %+.3f (dashed), %.0f%% positive. Preregistered.", S1_obs, 100 * S2_obs),
       x = NULL, y = "Kendall tau", colour = NULL) +
  theme(legend.position = "bottom", axis.text.x = element_blank())
ggsave(file.path(FIG, "P2_speccurve.png"), p_top, width = 10, height = 5, dpi = 150)

long <- do.call(rbind, lapply(c("series", "indicator", "detrend", "win", "spacing", "taper"),
  function(v) data.frame(dim = v, level = as.character(sc[[v]]), tau = sc$tau)))
p_dim <- ggplot(long, aes(reorder(level, tau, median), tau)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_boxplot(outlier.size = .3, fill = "grey90") +
  facet_wrap(~dim, scales = "free_x", nrow = 2) +
  labs(title = "Which analytic choice moves the answer?", x = NULL, y = "Kendall tau") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(FIG, "P2_by_dimension.png"), p_dim, width = 11, height = 7, dpi = 150)
