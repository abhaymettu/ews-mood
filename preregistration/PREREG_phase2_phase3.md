# Preregistration: Phases 2 and 3

**Project:** Early warning signals for critical transitions in mood
**Author:** Abhay Mettu
**Written:** 2026-07-26, before any Phase 2 or Phase 3 analysis was run.
**Commit:** the git commit that adds this file precedes every commit containing
Phase 2 or Phase 3 results. `git log --follow` on this file establishes the
ordering; that ordering is the only thing making this document a
preregistration rather than a description.

**Dataset:** Kossakowski, J. J., Groot, P. C., Haslbeck, J. M. B., Borsboom, D.,
& Wichers, M. (2017). Data from 'Critical Slowing Down as a Personalized Early
Warning Signal for Depression'. *Journal of Open Psychology Data*, 5: 1.
doi:10.5334/jopd.29. CC-BY 4.0.

---

## 0. What is NOT blind about this preregistration

This must come first, because pretending otherwise would be the exact practice
the project criticises.

**Phases 0 and 1 have already been run and their results are known to me.**
Specifically I already know:

1. On simulated data with a planted transition, no design in a 216-cell grid
   reached 80% power with a false-alarm rate below 10%.
2. A nominal Kendall p-value on overlapping rolling windows rejects 28-43% of
   the time under the null.
3. The variance indicator is confounded by where an item's mean sits relative
   to its response-category boundaries; autocorrelation indicators are not.
4. In this dataset, 13 of 16 pre-transition tests are significant with nominal
   p-values and 1 of 16 with a block-bootstrap null.
5. In this dataset, the variance signal is concentrated in the three items
   sitting on their floor (mean tau +0.62 vs -0.07 mid-scale), and a synthetic
   series with constant latent variance reproduces 75% of it. For 12 of 12
   items the observed trend is not beyond what rounding alone produces.

So this document does **not** preregister a blind test of whether early warning
signals are present. That question has effectively been answered for this
participant by Phase 1.

What it does preregister is the thing not yet known: **how much the conclusion
moves across the analyst degrees of freedom**, and what design follows from
that. The specification curve has not been computed. I do not know its median,
its spread, or what fraction of specifications is positive. Those are the
quantities registered below.

Where a prediction can still be made in advance, it is stated as a numbered
prediction with a falsification condition. Those are genuine: they are written
now and checked afterwards.

---

## 1. Phase 2 — specification curve

### 1.1 Question

Across the full space of defensible analytic choices, what is the distribution
of the early-warning trend statistic, and does the published conclusion depend
on a particular combination?

### 1.2 Analytic degrees of freedom to be crossed

Every combination of the following is run. No combination is dropped after
seeing its result.

| dimension | levels |
|---|---|
| rolling window width | 20, 30, 45, 60 days |
| detrending | none; linear; Gaussian kernel bw 7, 14, 28 days |
| series analysed | positive-affect composite; negative-affect composite; floored-item composite (suspic, doubt, irritat); mid-scale-item composite |
| overnight gap / spacing | native irregular; within-day pairs only (gap <= 3h); linear interpolation to a regular 3-hour grid |
| medication taper | ignored; included as a covariate (residualised); modelled as a structural break (phase-wise centring) |
| transition definition | changepoint on weekly SCL-90 (day 127); first sustained crossing of baseline mean + 2 SD (day 141) |
| indicator | variance; AC1 naive; AC1 within-day; AC1 continuous-time OU |

4 x 5 x 4 x 3 x 3 x 2 x 4 = **5,760 specifications.**

Per-item results (12 items) are reported separately as a descriptive
supplement, not folded into the curve, to avoid the item dimension dominating
it by sheer count.

### 1.3 Inference

Nominal p-values are not used anywhere, per DECISIONS.md D9.

**Per specification:** Kendall's tau of the rolling indicator against window
time, over the pre-transition window.

**Joint test (the registered inference):** a moving-block bootstrap applied to
the whole curve. For each of 100 surrogate datasets, blocks of the observed
series are resampled with replacement, original observation times reassigned,
and **the entire 5,760-specification curve recomputed**. This yields a null
distribution for curve-level summaries. Two summaries are registered:

- **S1** the median tau across all specifications
- **S2** the proportion of specifications with tau > 0

Block length is one day of beeps (median 6 observations); a 3-day block length
is run as a sensitivity analysis and both are reported.

### 1.4 Registered predictions

- **P1.** The median tau across the curve will be positive but its bootstrap
  p-value will exceed .05.
  *Falsified if* the joint test on S1 gives p <= .05.
- **P2.** Specifications using the variance indicator on the floored-item
  composite will occupy the top decile of the tau distribution.
  *Falsified if* fewer than half the top decile are variance-on-floored-items.
- **P3.** Detrending bandwidth will move tau more than window width does,
  measured as the range of median tau across levels of each.
  *Falsified if* window width produces the larger range.
- **P4.** No dimension will flip the sign of the median tau on its own except
  the series analysed.
  *Falsified if* any other dimension produces medians of opposite sign across
  its levels.

### 1.5 Context confound

Emotions shift with day of week and with reported events, which can create
false alarms and mask real ones. Registered test: recompute the curve after
residualising each series on (a) day of week, (b) the event-pleasantness and
event-importance items, (c) both. Report the shift in median tau and in the
timing of the maximum rolling indicator.

- **P5.** Adjusting for context will change the median tau by less than 0.10.
  *Falsified if* the change exceeds 0.10 in absolute value.

### 1.6 Stopping and reporting rules

- Every specification run is reported. There is no filtering step.
- Specifications that fail to compute (too few windows, non-convergent MLE) are
  reported as a count with reasons, not silently dropped.
- If the finding survives only a narrow region of the space, that region is
  named explicitly and reported as the headline.

---

## 2. Phase 3 — design recommendations

### 2.1 Question

Given a real ESM protocol, how much data is needed, and how must it be
processed, for critical slowing down to be a usable warning signal rather than
an artefact of analytic choices?

### 2.2 Inputs

Phase 0's power grid (216 cells, simulated, ground truth known) and Phase 2's
specification curve (real data, ground truth unknown). Neither alone answers
it: Phase 0 gives power without realism, Phase 2 gives realism without truth.

### 2.3 Deliverable

A design table stating, for each combination of study length, prompts per day,
tolerable missingness and response-scale resolution:

- detection power at a calibrated 5% level
- false alarm rate on a deteriorating-but-not-transitioning person
- **the false alarm rate a clinician would face**, expressed as expected false
  alarms per 100 monitored patients per year, which is the number that decides
  whether a monitoring system is deployable

### 2.4 Registered decision rule

A configuration is called **usable** only if it achieves power >= 0.80 *and*
false alarm rate <= 0.10. Both bars are conventions; they are fixed here so
they cannot be relaxed later to manufacture a positive recommendation.

- **P6.** No configuration within the ranges tested in Phase 0 will be usable
  under that rule, and the recommendation will therefore be conditional rather
  than a green light.
  *Falsified if* any configuration clears both bars.

*(P6 is a near-certain prediction: Phase 0 already established it for the
simulated grid. It is registered anyway so that Phase 3 cannot quietly relax
the rule when extending to the response-resolution dimension, which has not
been crossed with the power grid.)*

### 2.5 What Phase 3 will not claim

- It will not claim these numbers generalise beyond the model in `R/model.R`.
- It will not claim a recommendation validated on N=1 is validated.
- It will not present a design as usable on power alone.

---

## 3. Deviations

Any departure from this document is recorded in `preregistration/DEVIATIONS.md`
with the reason, and referenced from the final report. An empty deviations file
means none occurred.

## 4. OSF

This document is timestamped by its git commit. Posting it to OSF as a formal
registration is a separate, outward-facing step and is the author's to take; it
has deliberately not been done automatically.
