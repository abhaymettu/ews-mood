## sim/03_degrade.R -- Phase 0 step 3: what each layer of ESM realism does to
## the signal, one layer at a time, on a single path.
##
## This is diagnostic, not inferential. The point is to see WHERE the signal
## dies before the power analysis tells us HOW MUCH. One path, one seed, so the
## comparison is paired: every panel is the same underlying latent trajectory.

source("R/setup.R")
set.seed(SEED + 7)   # same path as sim/01_ideal.R

path <- simulate_latent(SIM_DAYS, "transition")
t_tr <- find_transition(path)
D    <- 180                                  # days of history before the event
cat(sprintf("transition at day %.1f; using the %d days before it\n", t_tr / 24, D))

## Layer stack. Each row adds exactly one thing to the row above it.
layers <- list(
  list(id = "1. continuous latent",   ppd = 48, sm = 0,   K = 0, mech = "none", rate = 0),
  list(id = "2. + 10 prompts/day",    ppd = 10, sm = 0,   K = 0, mech = "none", rate = 0),
  list(id = "3. + measurement noise", ppd = 10, sm = 0.4, K = 0, mech = "none", rate = 0),
  list(id = "4. + 7-point Likert",    ppd = 10, sm = 0.4, K = 7, mech = "none", rate = 0),
  list(id = "5. + 30% MCAR",          ppd = 10, sm = 0.4, K = 7, mech = "mcar", rate = 0.30),
  list(id = "6. + 30% MNAR instead",  ppd = 10, sm = 0.4, K = 7, mech = "mnar", rate = 0.30)
)

INDS <- c("variance", "ac1_naive", "ac1_withinday", "ac1_ou")
out <- list(); series <- list()

for (L in layers) {
  set.seed(SEED + 99)                        # same prompt schedule / noise draws
  o <- observe(path, end_hour = t_tr, window_days = D, prompts_per_day = L$ppd,
               sigma_m = L$sm, K = L$K, mechanism = L$mech, miss_rate = L$rate)
  pre <- detrend(o$y, o$t, "gaussian", bw_days = 14)
  for (k in INDS) {
    r <- rolling(pre, o$t, INDICATORS[[k]], win_days = 45, step_days = 2)
    if (is.null(r)) next
    r$indicator <- k; r$layer <- L$id
    out[[length(out) + 1]] <- r
  }
  o$layer <- L$id; o$n <- nrow(o)
  series[[length(series) + 1]] <- o
}
out    <- do.call(rbind, out)
series <- do.call(rbind, series)

nobs <- tapply(series$n, series$layer, function(z) z[1])
cat("\nobservations retained per layer:\n"); print(nobs)

## --- Kendall tau by layer, ACROSS PATHS --------------------------------------
## Deliberately not reported from the single path above. The between-path SD of
## tau under this design is ~0.33, which is larger than most of the layer
## effects we are trying to see: on one trajectory the layer ordering reverses
## routinely. This is a small, self-inflicted instance of exactly the Bos & De
## Jonge (2014) problem -- a single realisation cannot establish a pattern -- and
## it is worth hitting once in simulation, where the ground truth is known,
## before meeting it in real N = 1 data.
NPATH <- 150
tau_rep <- do.call(rbind, parallel::mclapply(1:NPATH, function(i) {
  set.seed(SEED + i)
  p <- simulate_latent(SIM_DAYS, "transition"); tt <- find_transition(p)
  if (is.na(tt) || tt < D * 24) return(NULL)     # need D days of history
  do.call(rbind, lapply(layers, function(L) {
    set.seed(SEED + 5000 + i)                    # paired across layers
    o <- observe(p, tt, D, L$ppd, sigma_m = L$sm, K = L$K,
                 mechanism = L$mech, miss_rate = L$rate)
    if (is.null(o)) return(NULL)
    pre <- detrend(o$y, o$t, "gaussian", bw_days = 14)
    z <- vapply(INDS, function(k) ews_tau(NULL, o$t, k, 45, 2, pre = pre), numeric(1))
    data.frame(rep = i, layer = L$id, t(z))
  }))
}, mc.cores = max(1L, parallel::detectCores() - 1L)))

agg <- function(f) {
  a <- aggregate(tau_rep[INDS], list(layer = tau_rep$layer), f)
  a[order(a$layer), ]
}
mn <- agg(function(z) mean(z, na.rm = TRUE)); sdv <- agg(function(z) sd(z, na.rm = TRUE))
cat(sprintf("\n--- mean Kendall tau over %d paths (SD across paths in brackets) ---\n",
            length(unique(tau_rep$rep))))
show <- mn
for (k in INDS) show[[k]] <- sprintf("%.3f (%.2f)", mn[[k]], sdv[[k]])
print(show, row.names = FALSE)
write.csv(merge(mn, sdv, by = "layer", suffixes = c("_mean", "_sd")),
          file.path(TAB, "03_degradation_tau.csv"), row.names = FALSE)

## --- the measurement-error attenuation, checked against theory ---------------
## An AR(1) observed with additive white noise has
##   plim(rho_hat) = rho * s2 / (s2 + sigma_m^2)
## so the naive AC1 is attenuated by the reliability ratio. As the system nears
## the fold, s2 grows, the ratio improves, and the naive AC1 rises PARTLY for
## that reason rather than because the recovery rate slowed. Here we measure
## both pieces on the same windows.
set.seed(SEED + 99)
o_clean <- observe(path, t_tr, D, 10, sigma_m = 0,   K = 0)
set.seed(SEED + 99)
o_noisy <- observe(path, t_tr, D, 10, sigma_m = 0.4, K = 0)
rc <- rolling(detrend(o_clean$y, o_clean$t, "gaussian", 14), o_clean$t, ind_ac1_naive, 45, 2)
rn <- rolling(detrend(o_noisy$y, o_noisy$t, "gaussian", 14), o_noisy$t, ind_ac1_naive, 45, 2)
vc <- rolling(detrend(o_clean$y, o_clean$t, "gaussian", 14), o_clean$t, ind_var, 45, 2)
att <- data.frame(day = rc$end / 24, clean = rc$val, noisy = rn$val,
                  s2 = vc$val)
att$predicted <- att$clean * att$s2 / (att$s2 + 0.4^2)
a <- na.omit(att)
a$R <- a$s2 / (a$s2 + 0.4^2)                       # reliability ratio
cat(sprintf("\nattenuation check: mean |observed noisy AC1 - theory| = %.4f (SD of noisy AC1 = %.3f)\n",
            mean(abs(a$noisy - a$predicted)), sd(a$noisy)))
## Exact decomposition of the observed rise: d(rho*R) = Rbar*d(rho) + rhobar*d(R)
d_obs <- diff(range_ends <- c(a$noisy[1], a$noisy[nrow(a)]))
d_rho <- a$clean[nrow(a)] - a$clean[1]
d_R   <- a$R[nrow(a)]   - a$R[1]
Rbar  <- mean(c(a$R[1], a$R[nrow(a)])); rhobar <- mean(c(a$clean[1], a$clean[nrow(a)]))
cat(sprintf("reliability ratio rises %.2f -> %.2f across the window\n", a$R[1], a$R[nrow(a)]))
cat(sprintf("observed AC1 rise %+.4f = %+.4f (true slowing) %+.4f (reliability improving)\n",
            d_obs, Rbar * d_rho, rhobar * d_R))
cat(sprintf("-> %.0f%% of the apparent 'rising autocorrelation' is measurement-error\n",
            100 * (rhobar * d_R) / (Rbar * d_rho + rhobar * d_R)))
cat("   attenuation shrinking, not the recovery rate slowing.\n")
write.csv(att, file.path(TAB, "03_attenuation.csv"), row.names = FALSE)

## --- figures -----------------------------------------------------------------
lab_i <- c(variance = "Variance", ac1_naive = "AC1 naive",
           ac1_withinday = "AC1 within-day", ac1_ou = "AC1 OU (irregular-aware)")
out$indicator <- factor(out$indicator, names(lab_i), lab_i)
out$day <- out$end / 24

p1 <- ggplot(out, aes(day, val, colour = layer)) +
  geom_line(linewidth = .45) +
  facet_wrap(~indicator, scales = "free_y") +
  scale_colour_viridis_d(end = .9) +
  labs(title = "Each layer of ESM realism, applied to the same latent path",
       subtitle = sprintf("45-day rolling windows, %d days before the transition", D),
       x = "day", y = NULL, colour = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "03_degradation_indicators.png"), p1, width = 10, height = 6.5, dpi = 150)

p2 <- ggplot(series, aes(t / 24, y)) +
  geom_point(size = .25, alpha = .5) +
  facet_wrap(~layer, ncol = 1, scales = "free_y") +
  labs(title = "The observed series at each layer", x = "day", y = "reported mood")
ggsave(file.path(FIG, "03_degradation_series.png"), p2, width = 9, height = 10, dpi = 150)

p3 <- ggplot(att, aes(day)) +
  geom_line(aes(y = clean, colour = "clean AC1"), linewidth = .5) +
  geom_line(aes(y = noisy, colour = "observed with noise"), linewidth = .5) +
  geom_line(aes(y = predicted, colour = "theory: rho * s2/(s2+sm2)"), linetype = 2) +
  scale_colour_manual(values = c("clean AC1" = "grey30", "observed with noise" = "steelblue4",
                                 "theory: rho * s2/(s2+sm2)" = "firebrick")) +
  labs(title = "Measurement error attenuates AC1, and the attenuation shrinks as the fold nears",
       subtitle = "part of the observed 'rising autocorrelation' is the reliability ratio improving",
       x = "day", y = "lag-1 autocorrelation", colour = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "03_attenuation.png"), p3, width = 8, height = 4.5, dpi = 150)

saveRDS(list(tau_mean = mn, tau_sd = sdv, tau_rep = tau_rep, att = att, nobs = nobs),
        file.path(TAB, "03_degrade.rds"))
