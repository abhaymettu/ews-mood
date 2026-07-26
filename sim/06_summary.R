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

## --- 1. the honest headline --------------------------------------------------
cat("=== 1. Best cell in the entire grid ===\n")
b <- pw[which.max(pw$power), ]
cat(sprintf("power %.2f at %d days x %d prompts/day, %s missingness, measurement noise %.1f,\n",
            b$power, b$days, b$ppd, b$miss, b$sigma_m))
cat(sprintf("indicator %s, window %.0f%% of series, %d observations.\n",
            b$indicator, 100 * b$win_frac, b$n_obs))
cat(sprintf("Its false alarm rate on a person who deteriorates without tipping: %.2f\n\n",
            b$fpr_drift))

## --- 2. minimum viable configurations ----------------------------------------
## For each measurement-noise / missingness combination, the cheapest design (by
## total prompts delivered) that clears 80% power.
pw$prompts_total <- round(pw$days * pw$ppd *
                            ifelse(pw$miss == "none", 1, 0.7))
viable <- subset(pw, power >= 0.80)
cat("=== 2. Cheapest design reaching 80% power, by noise and missingness ===\n")
if (!nrow(viable)) {
  cat("NONE. No configuration in the grid reaches 80% power.\n\n")
} else {
  mv <- do.call(rbind, lapply(split(viable, viable[c("sigma_m", "miss")], drop = TRUE),
                              function(z) z[which.min(z$prompts_total), ]))
  print(mv[order(mv$sigma_m, mv$miss),
           c("sigma_m", "miss", "days", "ppd", "indicator", "win_frac",
             "prompts_total", "power", "fpr_drift", "auc")],
        row.names = FALSE, digits = 3)
  cat("\n")
}

## --- 3. what the grid says you cannot buy ------------------------------------
cat("=== 3. Configurations where nothing works (max power over all indicators) ===\n")
cell <- do.call(rbind, lapply(split(pw, pw[c("days", "ppd", "miss", "sigma_m")], drop = TRUE),
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
for (v in c("days", "ppd", "miss", "sigma_m", "indicator", "win_frac")) {
  a <- aggregate(pw[c("power", "fpr_drift", "auc")], list(pw[[v]]), mean, na.rm = TRUE)
  names(a)[1] <- v
  print(a, row.names = FALSE, digits = 3); cat("\n")
}

## --- 5. days vs prompts/day: which buys more? --------------------------------
## Both cost the participant. Holding total prompts roughly constant, is it
## better to sample densely for a short time or sparsely for a long time?
cat("=== 5. Same budget, spent two ways (no missingness, all indicators pooled) ===\n")
nm <- subset(pw, miss == "none")
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
