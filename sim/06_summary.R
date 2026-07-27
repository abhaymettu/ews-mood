## sim/06_summary.R -- the Phase 0 deliverable: what is the minimum viable series?
##
## Reads the tables produced by 02, 04 and 05 and answers the question the power
## analysis was built for. Does no simulation of its own.
##
## "Viable" is defined here as power >= 0.80 at a 5% level calibrated against a
## no-change null. That is a conventional bar, not a clinical one; the false
## alarm column is what a clinician would actually care about and it is reported
## alongside every configuration rather than buried.

source("R/setup.R")

pw  <- read.csv(file.path(TAB, "04_power.csv"))
thr <- read.csv(file.path(TAB, "02_thresholds.csv"))

pw$resp <- ifelse(pw$K == 0, "continuous", "Likert-7")
pw$prompts_total <- round(pw$days * pw$ppd * ifelse(pw$miss == "none", 1, 0.7))

## --- 1. the honest headline --------------------------------------------------
## Power alone is not a headline. A detector that fires on everyone has power 1.
## The pair that matters is (power, false alarm rate on someone who deteriorates
## and does not become ill).
cat("=== 1. Highest-power cell, and what it costs in false alarms ===\n")
b <- pw[which.max(pw$power), ]
cat(sprintf("power %.2f at %d days x %d prompts/day, %s missingness, sigma_m %.1f, %s response\n",
            b$power, b$days, b$ppd, b$miss, b$sigma_m, b$resp))
cat(sprintf("indicator %s, %d observations.\n", b$indicator, b$n_obs))
cat(sprintf("FALSE ALARM RATE on a deteriorating-but-not-tipping person: %.2f\n\n", b$fpr_drift))

## --- 2. minimum viable configurations ----------------------------------------
## Viable = 80% power AND no more than 10% false alarms. Both bars are
## conventions, but dropping either one makes the answer meaningless.
BAR_POWER <- 0.80; BAR_FPR <- 0.10
cat(sprintf("=== 2. Cheapest design with power >= %.0f%% AND false alarms <= %.0f%% ===\n",
            100 * BAR_POWER, 100 * BAR_FPR))
viable <- subset(pw, power >= BAR_POWER & fpr_drift <= BAR_FPR)
if (!nrow(viable)) {
  cat("NONE. No configuration in the grid clears both bars.\n")
  cat("\nCells clearing the power bar alone, with the false alarm rate they carry:\n")
  vp <- subset(pw, power >= BAR_POWER)
  if (!nrow(vp)) cat("  (none clear even the power bar)\n") else
    print(vp[order(vp$fpr_drift),
             c("resp", "days", "ppd", "miss", "sigma_m", "indicator",
               "prompts_total", "power", "fpr_drift", "auc")],
          row.names = FALSE, digits = 3)
  cat("\nCells clearing the false-alarm bar alone, best power:\n")
  vf <- subset(pw, fpr_drift <= BAR_FPR)
  print(head(vf[order(-vf$power),
                c("resp", "days", "ppd", "miss", "sigma_m", "indicator",
                  "prompts_total", "power", "fpr_drift", "auc")], 8),
        row.names = FALSE, digits = 3)
  cat("\n")
} else {
  mv <- do.call(rbind, lapply(split(viable, viable[c("sigma_m", "miss", "K")], drop = TRUE),
                              function(z) z[which.min(z$prompts_total), ]))
  print(mv[order(mv$sigma_m, mv$miss),
           c("resp", "sigma_m", "miss", "days", "ppd", "indicator",
             "prompts_total", "power", "fpr_drift", "auc")],
        row.names = FALSE, digits = 3)
  cat("\n")
}

## --- 2b. the artefact, in the grid -------------------------------------------
## sim/07 showed the Likert rounding artefact inflates the variance indicator.
## Here is the same thing inside the full design grid: how much of each
## indicator's discrimination survives going continuous.
cat("=== 2b. How much of each indicator's performance is the response scale? ===\n")
rk <- aggregate(pw[c("power", "fpr_drift", "auc")], pw[c("indicator", "resp")],
                mean, na.rm = TRUE)
w <- reshape(rk, idvar = "indicator", timevar = "resp", direction = "wide")
w$auc_inflation <- w$`auc.Likert-7` - w$auc.continuous
print(w[order(-w$auc_inflation), ], row.names = FALSE, digits = 3)
cat("\nAC indicators are unmoved by discretisation. Variance is not.\n")
cat("Its lead over the best AC indicator on continuous data is the honest one.\n\n")

## --- 3. what the grid says you cannot buy ------------------------------------
cat("=== 3. Configurations where nothing works (max power over all indicators) ===\n")
cell <- do.call(rbind, lapply(split(pw, pw[c("days", "ppd", "miss", "sigma_m", "K")], drop = TRUE),
                              function(z) z[which.max(z$power), ]))
dead <- subset(cell, power < 0.30)
cat(sprintf("%d of %d design cells cannot reach even 30%% power with any indicator.\n",
            nrow(dead), nrow(cell)))
if (nrow(dead)) {
  cat("By measurement noise:\n"); print(table(dead$sigma_m))
  cat("By days:\n"); print(table(dead$days))
}

## --- 4. marginal cost of each design lever -----------------------------------
cat("\n=== 4. Marginal effect of each lever (mean over the rest of the grid) ===\n")
for (v in c("days", "ppd", "miss", "sigma_m", "resp", "indicator")) {
  a <- aggregate(pw[c("power", "fpr_drift", "auc")], list(pw[[v]]), mean, na.rm = TRUE)
  names(a)[1] <- v
  print(a, row.names = FALSE, digits = 3); cat("\n")
}

## --- 5. days vs prompts/day: which buys more? --------------------------------
## Both cost the participant. Holding total prompts roughly constant, is it
## better to sample densely for a short time or sparsely for a long time?
cat("=== 5. Same budget, spent two ways (no missingness, all indicators pooled) ===\n")
nm <- subset(pw, miss == "none" & resp == "continuous")
tradeoff <- aggregate(nm[c("power", "auc")], nm[c("days", "ppd")], max, na.rm = TRUE)
tradeoff$prompts <- tradeoff$days * tradeoff$ppd
print(tradeoff[order(tradeoff$prompts), ], row.names = FALSE, digits = 3)

## --- 6. the size of the analytic-p-value error -------------------------------
cat("\n=== 6. Cost of using a nominal Kendall p-value instead of a calibrated one ===\n")
cat(sprintf("nominal 5%% critical value:   tau = %.2f to %.2f\n",
            min(thr$thr_nominal), max(thr$thr_nominal)))
cat(sprintf("actual 5%% point of the null: tau = %.2f to %.2f\n",
            min(thr$thr_emp95), max(thr$thr_emp95)))
cat(sprintf("=> a nominally 5%% test actually rejects %.0f-%.0f%% of the time under the null.\n",
            100 * min(thr$fpr_if_nominal), 100 * max(thr$fpr_if_nominal)))

## Redirect with the shell if you want this on disk:
##   Rscript sim/06_summary.R | tee output/tables/06_summary.txt
