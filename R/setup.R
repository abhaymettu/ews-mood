## R/setup.R -- sourced by every script. Seed, paths, plot theme.

SEED <- 20260726L   # recorded in README; every script re-seeds from it

## Shared simulation geometry. Every script uses these so that figures, null
## calibration and the power grid all describe the same generative process.
##   HOLD_DAYS  c is held flat this long before the ramp begins, so that even a
##              180-day analysis window has history in front of it for EVERY
##              path rather than only for the late-tipping ones (DECISIONS.md D13)
##   SIM_DAYS   HOLD_DAYS + 300 days of ramp; the 300 sets the approach speed
##   REF_DAY    median tipping day under these settings; null observation
##              windows end here so they sit at the same calendar position as
##              the transition windows they are compared against
HOLD_DAYS <- 120
SIM_DAYS  <- 420
REF_DAY   <- 314

## Run everything from the project root. One assertion beats path-sniffing.
PROJ <- getwd()
stopifnot("run scripts from the project root" = dir.exists(file.path(PROJ, "R")))

FIG <- file.path(PROJ, "output", "figures")
TAB <- file.path(PROJ, "output", "tables")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

for (f in c("model.R", "observe.R", "detect.R"))
  source(file.path(PROJ, "R", f))

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")))

set.seed(SEED)
