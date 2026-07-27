## analysis/05_design.R -- Phase 3: design recommendations
##
## Preregistered in preregistration/PREREG_phase2_phase3.md section 2.
##
## The question: given a real ESM protocol, how much data do you need, and how
## must it be processed, for critical slowing down to be a usable warning signal
## rather than an artefact of analytic choices?
##
## Inputs. Phase 0's power grid (simulated, ground truth known, 216 cells) and
## Phase 2's specification curve (real data, ground truth unknown). Neither
## alone answers it: Phase 0 has power without realism, Phase 2 has realism
## without truth.
##
## THE STEP THAT MATTERS. Power and false-alarm rate are per-assessment
## quantities. A clinician does not experience per-assessment quantities; they
## experience alarms per patient per year, and they act on the probability that
## an alarm is real. That depends on the base rate of transitions, which no
## amount of statistical power can change. Translating into those units is the
## whole point of this phase, and it is where an econometrics reflex --
## always ask what the decision-maker actually faces -- earns its keep.
##
## Registered usability rule, fixed in advance: power >= 0.80 AND
## false alarm rate <= 0.10.

source("R/setup.R")

pw <- read.csv(file.path(TAB, "04_power.csv"))
pw$resp <- ifelse(pw$K == 0, "continuous", "Likert-7")

## --- deployment translation --------------------------------------------------
## A D-day monitoring window produces one assessment. Monitoring a patient for a
## year with non-overlapping windows gives 365/D assessments. Non-overlapping is
## the conservative choice: sliding windows produce more alarms, but correlated
## ones, and modelling that correlation would need assumptions we cannot support.
pw$assess_per_year <- 365 / pw$days
pw$false_alarms_per_100py <- 100 * pw$fpr_drift * pw$assess_per_year

## Base rate of a genuine transition per patient-year. Swept, because it moves
## the answer more than anything the analyst controls, and because the right
## value depends entirely on who is being monitored:
##   0.05  roughly general-population annual incidence of a depressive episode
##   0.20  elevated-risk / subclinical
##   0.40  remitted MDD, which is the Kossakowski participant's situation
##         (a patient tapering antidepressants has a high relapse risk)
BASE <- c(0.05, 0.20, 0.40)

ppv_at <- function(power, fpr, days, base) {
  n_assess <- 365 / days
  ## a transitioning patient gets one chance to be caught correctly per year
  tp <- base * power
  ## a non-transitioning patient can false-alarm at each assessment
  fp <- (1 - base) * (1 - (1 - fpr)^n_assess)
  tp / (tp + fp)
}
for (b in BASE) pw[[sprintf("ppv_%02d", round(100 * b))]] <-
  ppv_at(pw$power, pw$fpr_drift, pw$days, b)

## --- registered usability rule ----------------------------------------------
BAR_POWER <- 0.80; BAR_FPR <- 0.10
pw$usable <- pw$power >= BAR_POWER & pw$fpr_drift <= BAR_FPR

cat("=== Registered rule: power >= 0.80 AND false alarm rate <= 0.10 ===\n")
cat(sprintf("usable configurations: %d of %d\n\n", sum(pw$usable), nrow(pw)))
cat(sprintf("P6 (no configuration is usable): %s\n\n",
            ifelse(sum(pw$usable) == 0, "SUPPORTED", "FALSIFIED")))

## --- the best you can do, and what it would mean in a clinic -----------------
cat("=== Best-power configurations, translated into deployment units ===\n")
top <- head(pw[order(-pw$power), ], 8)
print(top[c("resp", "days", "ppd", "miss", "sigma_m", "indicator", "power",
            "fpr_drift", "false_alarms_per_100py", "ppv_05", "ppv_20", "ppv_40")],
      row.names = FALSE, digits = 3)

cat("\nRead the last three columns as: of every 100 alarms raised, this many are real.\n")
cat("Base rates: 5% general population, 20% elevated risk, 40% remitted MDD tapering\n")
cat("medication (the situation of the participant in the reproduction dataset).\n")

## --- what would it take? -----------------------------------------------------
## Invert the question: holding the design fixed at the best cell, what false
## alarm rate would be needed for a majority of alarms to be real?
b <- pw[which.max(pw$power), ]
need <- vapply(BASE, function(base) {
  f <- uniroot(function(fp) ppv_at(b$power, fp, b$days, base) - 0.5,
               c(1e-6, 0.99))$root
  f
}, numeric(1))
cat(sprintf("\n=== What the best cell (%d days, %d/day, %s) would need ===\n",
            b$days, b$ppd, b$resp))
cat(sprintf("achieved false alarm rate: %.3f\n", b$fpr_drift))
for (i in seq_along(BASE))
  cat(sprintf("  to make >50%% of alarms real at a %.0f%% base rate, it needs <= %.4f  (%.0fx better)\n",
              100 * BASE[i], need[i], b$fpr_drift / need[i]))

## --- the design table --------------------------------------------------------
## Best indicator per design cell, restricted to the continuous-response arm,
## because Phase 1 showed the variance indicator on coarse scales is not
## measuring what it claims to (analysis/03_artefact_test.R).
cont <- subset(pw, resp == "continuous")
cell <- do.call(rbind, lapply(split(cont, cont[c("days", "ppd", "miss", "sigma_m")],
                                    drop = TRUE), function(z) z[which.max(z$power), ]))
cell <- cell[order(-cell$power), ]
write.csv(cell, file.path(TAB, "P3_design_table.csv"), row.names = FALSE)

cat("\n=== Design table (continuous responses, best indicator per cell) ===\n")
print(head(cell[c("days", "ppd", "miss", "sigma_m", "indicator", "power",
                  "fpr_drift", "false_alarms_per_100py", "ppv_40")], 12),
      row.names = FALSE, digits = 3)

## --- marginal guidance -------------------------------------------------------
cat("\n=== What each design lever buys (continuous arm) ===\n")
for (v in c("days", "ppd", "miss", "sigma_m")) {
  a <- aggregate(cont[c("power", "fpr_drift", "false_alarms_per_100py", "ppv_40")],
                 list(cont[[v]]), mean, na.rm = TRUE)
  names(a)[1] <- v
  print(a, row.names = FALSE, digits = 3); cat("\n")
}

## --- duration vs density, the one clear recommendation -----------------------
nm <- subset(cont, miss == "none")
td <- aggregate(nm[c("power", "auc")], nm[c("days", "ppd")], max, na.rm = TRUE)
td$prompts <- td$days * td$ppd
td$per_prompt <- td$power / td$prompts * 1000
cat("=== Duration vs density at equal budget ===\n")
print(td[order(-td$per_prompt), ], row.names = FALSE, digits = 3)

saveRDS(list(pw = pw, cell = cell, need = need, base = BASE),
        file.path(TAB, "P3_design.rds"))

## --- figures -----------------------------------------------------------------
long <- do.call(rbind, lapply(seq_along(BASE), function(i)
  data.frame(days = cont$days, ppd = cont$ppd, power = cont$power,
             ppv = cont[[sprintf("ppv_%02d", round(100 * BASE[i]))]],
             base = sprintf("base rate %.0f%%", 100 * BASE[i]))))
p1 <- ggplot(long, aes(power, ppv, colour = factor(days))) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_point(alpha = .7, size = 1) +
  facet_wrap(~base) +
  scale_colour_viridis_d(end = .85, name = "days") +
  labs(title = "Detection power buys much less than it looks like",
       subtitle = "y = share of raised alarms that are real. Dashed line = half. Simulated grid, continuous responses.",
       x = "power", y = "positive predictive value")
ggsave(file.path(FIG, "P3_ppv.png"), p1, width = 10, height = 4.5, dpi = 150)

p2 <- ggplot(td, aes(days, ppd, fill = power)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", power)), size = 3) +
  scale_fill_viridis_c(option = "magma", end = .9) +
  labs(title = "Power by study length and sampling density",
       subtitle = "continuous responses, no missingness, best indicator",
       x = "days of monitoring", y = "prompts per day")
ggsave(file.path(FIG, "P3_duration_density.png"), p2, width = 7, height = 4.5, dpi = 150)
