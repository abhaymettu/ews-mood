## analysis/01_clean.R -- clean the Kossakowski et al. (2017) ESM data
##
## Every choice here is contestable and every one is logged to DECISIONS.md
## (D15-D21). Phase 2 sweeps the ones that can be swept. Nothing in this script
## looks at any early-warning statistic: the cleaning is fixed before the
## analysis touches it.
##
## Data: Kossakowski et al. (2017), JOPD 5:1, doi 10.5334/jopd.29, CC-BY.

source("R/setup.R")

raw <- read.csv(file.path(PROJ, "data-raw", "ESMdata", "ESMdata.csv"),
                stringsAsFactors = FALSE)
log_lines <- c()
note <- function(...) {
  s <- sprintf(...); log_lines <<- c(log_lines, s); cat(s, "\n")
}

## --- D15: build a real time axis --------------------------------------------
## `dayno` in the raw file is day-of-YEAR, so it runs 226..366 then wraps to
## 1..98 when the study crosses New Year 2012/13. Using it as a time index
## silently folds the second half of the study back on top of the first. The
## time axis is rebuilt from `date` + `beeptime` instead.
d <- raw
d$date_p <- as.Date(d$date, format = "%d/%m/%y")
stopifnot("unparsed dates" = !any(is.na(d$date_p)))
hm <- do.call(rbind, strsplit(d$beeptime, ":"))
d$t_hours <- as.numeric(d$date_p - min(d$date_p)) * 24 +
             as.numeric(hm[, 1]) + as.numeric(hm[, 2]) / 60
d$day <- as.numeric(d$date_p - min(d$date_p)) + 1
d <- d[order(d$t_hours), ]
note("D15 time axis rebuilt from date+beeptime: %d days, %.0f hours span",
     max(d$day), diff(range(d$t_hours)))
note("     (raw dayno is day-of-year and wraps at New Year; not used)")

## --- D16: the mood items are on two different scales -------------------------
## Eight items run 1..7; four negative-valence items (down, lonely, anxious,
## guilty) run -3..+3. Both are 7-point, differently anchored.
##
## They are NOT rescaled to a common range. Rescaling is a linear map, so it
## leaves autocorrelation untouched and multiplies variance by a constant --
## which changes nothing for a rank-based trend statistic like Kendall's tau.
## What it would destroy is the relationship between an item's mean and its
## rounding boundaries, and Phase 0 (PHASE0.md section 4) showed that
## relationship drives a large artefact in the variance indicator. Keeping the
## items on their native scales keeps that artefact visible rather than
## smearing it.
POS <- c("mood_relaxed", "mood_satisfi", "mood_enthus", "mood_cheerf", "mood_strong")
NEG <- c("mood_down", "mood_lonely", "mood_anxious", "mood_guilty",
         "mood_irritat", "mood_suspic", "mood_doubt")
MOOD <- c(POS, NEG)

rng <- t(sapply(MOOD, function(m) c(min = min(d[[m]], na.rm = TRUE),
                                    max = max(d[[m]], na.rm = TRUE),
                                    mean = mean(d[[m]], na.rm = TRUE))))
## Distance from the item mean to its nearest response boundary, in categories.
## Phase 0 predicts items sitting near a floor or ceiling behave differently
## from items sitting mid-scale. This is the variable that prediction is about.
scale_lo <- ifelse(MOOD %in% c("mood_down", "mood_lonely", "mood_anxious", "mood_guilty"), -3, 1)
scale_hi <- scale_lo + 6
rng <- data.frame(item = MOOD, rng,
                  headroom_lo = rng[, "mean"] - scale_lo,
                  headroom_hi = scale_hi - rng[, "mean"])
note("D16 %d mood items on two anchorings (1..7 and -3..+3); left on native scales",
     length(MOOD))
note("     items within 1.5 categories of a floor/ceiling: %s",
     paste(rng$item[pmin(rng$headroom_lo, rng$headroom_hi) < 1.5], collapse = ", "))

## --- D17: composites ---------------------------------------------------------
## Both a positive-affect and a negative-affect composite, plus every item kept
## individually. The published analyses use composites; Phase 2 sweeps
## composite-vs-item because it is a live analyst degree of freedom.
d$pa <- rowMeans(d[POS], na.rm = TRUE)
d$na_ <- rowMeans(d[NEG], na.rm = TRUE)
note("D17 composites: pa = mean(%d positive items), na_ = mean(%d negative items)",
     length(POS), length(NEG))

## --- D18: aborted responses --------------------------------------------------
## resp_abort flags beeps the participant started and did not finish. 7 of 1476.
## Dropped: a partially completed questionnaire is not a measurement of mood at
## that moment, and 7 rows cannot matter either way. Recorded so it is not
## silent.
n0 <- nrow(d)
d <- d[is.na(d$resp_abort) | d$resp_abort == 0, ]
note("D18 dropped %d aborted responses (%.1f%%)", n0 - nrow(d), 100 * (n0 - nrow(d)) / n0)

## --- D19: missingness --------------------------------------------------------
## Item-level missingness is negligible (<0.3%). The real missingness is beeps
## that never happened: the protocol was 10/day, the median observed is 6.
## Rows are NOT imputed here -- the detectors in R/detect.R handle irregular
## spacing natively, and imputation is one of the analyst degrees of freedom
## Phase 2 sweeps. Rows with a missing composite are dropped.
n0 <- nrow(d)
d <- d[!is.na(d$pa) & !is.na(d$na_), ]
bpd <- table(d$day)
note("D19 dropped %d rows with a missing composite; %d beeps over %d days",
     n0 - nrow(d), nrow(d), length(bpd))
note("     beeps/day: median %d, IQR %d-%d, range %d-%d (protocol was 10)",
     median(bpd), quantile(bpd, .25), quantile(bpd, .75), min(bpd), max(bpd))
note("     realised compliance %.0f%% of a 10/day protocol", 100 * nrow(d) / (10 * max(d$day)))

## --- D20: study phases and the medication taper ------------------------------
## 1 baseline | 2 double-blind, no reduction | 3 double-blind dose reduction
## 4 post-assessment | 5 follow-up.  The taper is the experimental manipulation
## and a candidate structural break. Coded here, modelled in Phase 2 three ways
## (covariate / structural break / ignored) because that choice is contestable.
ph <- do.call(rbind, lapply(split(d$day, d$phase), function(z)
  data.frame(first = min(z), last = max(z), days = max(z) - min(z) + 1)))
ph$phase <- rownames(ph)
d$tapering <- as.integer(d$phase == 3)
d$post_taper <- as.integer(d$phase >= 4)
note("D20 phase spans (day): %s",
     paste(sprintf("p%s=%d-%d", ph$phase, ph$first, ph$last), collapse = "  "))

## --- D21: when is the transition? -------------------------------------------
## This is the most consequential decision in Phase 1 and the one the published
## framing is loosest about.
##
## The weekly SCL-90 depression score does NOT jump. It sits at 1.0-1.6 through
## phases 1-3 and rises through phase 4 to a sustained 1.8-2.5 in phase 5. That
## is a gradual escalation, not a visible discontinuity, and calling it a
## "critical transition" is an interpretation, not an observation.
##
## We locate it two ways and carry BOTH forward, because they disagree:
##   (a) least-squares changepoint on the weekly score (objective, no threshold)
##   (b) first sustained crossing of baseline mean + 2 SD, baseline = phases 1-3
scl <- unique(d[!is.na(d$dep), c("day", "phase", "dep")])
scl <- scl[order(scl$day), ]

ss <- vapply(2:(nrow(scl) - 1), function(k) {
  a <- scl$dep[1:k]; b <- scl$dep[(k + 1):nrow(scl)]
  sum((a - mean(a))^2) + sum((b - mean(b))^2)
}, numeric(1))
cp_day <- scl$day[(2:(nrow(scl) - 1))[which.min(ss)]]

base <- scl$dep[scl$phase %in% 1:3]
thr <- mean(base) + 2 * sd(base)
above <- scl$dep > thr
sust <- which(above & c(above[-1], FALSE))          # crossing that stays up
cross_day <- if (length(sust)) scl$day[sust[1]] else NA_integer_

note("D21 baseline (phases 1-3) SCL-90: mean %.2f SD %.2f, threshold %.2f", mean(base), sd(base), thr)
note("     (a) changepoint on weekly score: day %d", cp_day)
note("     (b) first sustained crossing of mean+2SD: day %d", cross_day)
note("     BOTH carried forward; Phase 2 sweeps which one is used")
note("     NOTE the score rises gradually. 'Transition' is an interpretation.")

TRANSITION <- list(changepoint = cp_day, crossing = cross_day)

## --- write -------------------------------------------------------------------
dir.create(file.path(PROJ, "data"), showWarnings = FALSE)
keep <- c("day", "t_hours", "date", "phase", "beepno", "tapering", "post_taper",
          MOOD, "pa", "na_", "dep")
clean <- d[keep]
write.csv(clean, file.path(PROJ, "data", "esm_clean.csv"), row.names = FALSE)
saveRDS(list(data = clean, transition = TRANSITION, mood_items = MOOD,
             pos = POS, neg = NEG, item_ranges = rng, phases = ph, scl = scl),
        file.path(PROJ, "data", "esm_clean.rds"))

writeLines(c("# Cleaning log -- analysis/01_clean.R",
             paste0("# generated ", format(Sys.time(), "%Y-%m-%d %H:%M")),
             "", log_lines),
           file.path(PROJ, "data", "cleaning_log.txt"))

cat("\n--- item scale positions (Phase 0 predicts these matter) ---\n")
print(rng[order(pmin(rng$headroom_lo, rng$headroom_hi)), ], row.names = FALSE, digits = 3)
cat(sprintf("\nwrote data/esm_clean.csv: %d rows x %d cols\n", nrow(clean), ncol(clean)))
