# Early warning signals for critical transitions in mood

Can rising autocorrelation and variance in intensive longitudinal mood data
warn of a transition into depression, in a single person, before it happens?

The hypothesis (van de Leemput et al., *PNAS* 2014; Wichers et al.,
*Psychother Psychosom* 2016) is that mood is a dynamical system with alternative
stable states, and that critical slowing down precedes a transition between
them. The published critique (Bos & De Jonge, *PNAS* 2014, conceded by the
original authors) is that between-subject evidence does not establish a
within-person phenomenon. Whether intraindividual early warning signals
anticipate intraindividual transitions is open, and it is a time-series problem
about short, noisy, unevenly-spaced series.

**Status: Phase 0 complete. No real data has been downloaded.**

Headline: on simulated data where the transition is planted and known, **no design
in a 216-cell grid reached 80% power with a false-alarm rate under 10%**. Three
separate mechanisms bias the standard analysis toward finding critical slowing
down whether or not it is present. Full write-up in [PHASE0.md](PHASE0.md);
every contestable choice in [DECISIONS.md](DECISIONS.md).

## Phases

| Phase | Question | Status |
|---|---|---|
| 0 | Does the detector work when the ground truth is known? | done, see PHASE0.md |
| 1 | Can the published finding be reproduced? (Kossakowski et al. 2017) | not started |
| 2 | How much does the conclusion move across analyst degrees of freedom? | not started, **preregister first** |
| 3 | How much data does a real ESM study need for this to be usable? | not started, **preregister first** |

Phase 0 exists so that no detector is ever tuned against the real data. Nothing
in `/data-raw` until Phase 0 is reviewed and signed off.

## Reproducing

```
Rscript sim/01_ideal.R              #  ~20 s   gate: CSD on clean data
Rscript sim/02_null_calibration.R   #  ~2 min  calibrated thresholds
Rscript sim/03_degrade.R            #  ~1 min  degradation + attenuation
Rscript sim/04_power.R              # ~11 min  the grid (cached; REFRESH=TRUE to redo)
Rscript sim/05_measurement_error.R  # ~28 min  state-space comparison
Rscript sim/07_likert_artifact.R    #  ~2 min  response-scale artefact
Rscript sim/06_summary.R            #   ~1 s   the design-recommendation tables
```

Run from the project root; the scripts assert this. Global seed `20260726`, set
in `R/setup.R` and re-set per replicate so results are invariant to core count.

## Layout

```
R/          model.R    stochastic double-well + control parameter
            observe.R  ESM sampling, measurement noise, Likert, missingness
            detect.R   detrending, EWS indicators, rolling windows, trend test
            setup.R    seed, paths, sourced by everything
sim/        01_ideal.R             gate: does CSD appear on clean data?
            02_null_calibration.R  what does "significant" mean here?
            03_degrade.R           what each layer of realism costs
            04_power.R             the power analysis
            05_measurement_error.R does modelling measurement error help? (no)
            06_summary.R           reads the tables, answers the design question
            07_likert_artifact.R   the response scale can reverse the signal
output/     figures/, tables/  (script-generated, safe to delete)
PHASE0.md      findings
DECISIONS.md   every contestable choice, what else was available, why this one
```

`data-raw/`, `data/`, `analysis/`, `preregistration/`, `report/` are empty and
belong to Phases 1-3.

## The model

```
dx = (-x^3 + x + c(t)) dt + sigma_p dW
```

The normal form of a fold bifurcation. `x` is latent mood valence, `c` a slow
control parameter. Two stable branches exist for `|c| < 0.3849`; lowering `c`
past `-0.3849` destroys the euthymic branch and the system drops into the
depressed one.

Critical slowing down is *derived*, not assumed: the local recovery rate is
`lambda = -3x*^2 + 1`, and on the euthymic branch `x* -> 1/sqrt(3)` as the fold
is approached, so `lambda -> 0`, the AR(1) coefficient `exp(lambda * gap) -> 1`
and the stationary variance `sigma_p^2/(2|lambda|) -> inf`. That derivation is
also the yardstick the detectors are checked against.

One model time unit is calibrated to 5 hours, so that far from the bifurcation
the lag-1 autocorrelation at 90-minute spacing is ~0.4, a typical ESM
"emotional inertia" value. Every statement about days of data depends on this
anchor (DECISIONS.md D3).

## Indicators

| | what it is | why it is here |
|---|---|---|
| `variance` | rolling variance of the detrended series | the standard indicator; robust to additive white noise but **badly confounded by Likert rounding** — see PHASE0.md section 4 |
| `ac1_naive` | `cor(y[-n], y[-1])` | what the applied literature computes; treats consecutive *observations* as lag-1 regardless of elapsed time |
| `ac1_withinday` | same, restricted to gaps <= 3 h | cheap fix for the overnight gap; costs sample size |
| `ac1_ou` | continuous-time OU fitted by exact ML on irregular gaps | estimates `lambda` directly, so irregular spacing is data rather than a nuisance |
| `ac1_ou_me` | OU + explicit measurement-error variance, Kalman filter | the only one that separates true slowing from shrinking attenuation — and it still loses, see PHASE0.md section 5 |

Hand-rolled rather than taken from `earlywarnings`, which requires an evenly
spaced `ts` object (DECISIONS.md D8).

## Environment

R 4.5.1. Phase 0 uses base R plus `ggplot2` only. That is deliberate: every
statistical operation here is `cor`, `optim`, `lm` and `dnorm`, and the point of
Phase 0 is to know exactly what the detector does.

`earlywarnings`, `mlVAR`, `graphicalVAR`, `qgraph`, `psychonetrics`, `brms`,
`osfr` and `quarto` are **not** installed yet. They are added in the phase that
first needs them, and `renv.lock` is snapshotted at that point, rather than
pinning a large stack that Phase 0 does not exercise.

## Relation to prior work

**Read this before the results.** A literature check was run after the analyses
and before publication, and it narrowed what this project can claim.

- The floor-effect objection tested here was raised against this dataset and
  **answered in the source paper** (Wichers et al. 2016). Their rebuttal holds
  for the autocorrelation indicators and fails for variance — a refinement, not
  a refutation.
- Helmich et al. (2024, *Nature Reviews Psychology*) **already conclude** that
  early warning signals do not reliably predict symptom change. This report
  converges with a published review; it is corroboration, not news.
- von Klipstein et al. (2023) showed floor effects explain a different ESM
  finding, so the general concern is established in the field.

What is added is narrow: a quantitative test where the published answer was
qualitative, the finding that the confound is specific to the variance
indicator, and the deployment arithmetic (PPV at realistic base rates). I could
not obtain full text of Helmich et al. (2024) — **treat the overlap with it as
unresolved.**

## Replication failed

The discretisation artefact was taken to a second dataset (Fisher et al. 2017,
40 participants, 840 item-series, 0-100 visual analog scales). Rounding those
same series to 7 categories **did not** create the pattern found in
Kossakowski: the near-bound split was -0.090 on VAS and -0.085 after
discretising (paired p = 0.59 across participants). The obvious rescue is ruled
out — Fisher's item means travel 3.5x *further* than the Kossakowski
participant's.

So the artefact is demonstrable in simulation and consistent with the
Kossakowski data, but **it is not a general property of coarse ESM scales.**
Whether the Kossakowski pattern reflects that participant's monotone
well-to-ill progression or is chance in one dataset is **not resolved here.**

## Scope

Independent work, not affiliated with or endorsed by the authors of any paper
discussed. What is criticised is a **method** — a trend statistic that is
anticonservative on overlapping windows, a measurement-error bias pointing
toward the hypothesis, and a discretisation artefact that can produce the
primary indicator's signal outright — not a research group. This analysis was
only possible because the original authors published item-level raw data under
CC-BY, which almost nobody does.

Nothing here shows critical slowing down in mood is absent. See
[report/report.qmd](report/report.qmd) sections on Errors, Limitations, and
Acknowledgement.

## Citation

The dataset used from Phase 1 onward is CC-BY and must be cited wherever it
appears:

> Kossakowski, J. J., Groot, P. C., Haslbeck, J. M. B., Borsboom, D., & Wichers,
> M. (2017). Data from 'Critical Slowing Down as a Personalized Early Warning
> Signal for Depression'. *Journal of Open Psychology Data*, 5: 1.
> https://doi.org/10.5334/jopd.29
