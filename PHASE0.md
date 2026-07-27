# Phase 0 findings: can we detect critical slowing down when we know it is there?

Simulation only. No real data has been downloaded. Every number below comes from
`sim/*.R` and regenerates from a clean clone; seed 20260726.

The question Phase 0 answers is deliberately narrow and deliberately easy:
**we planted a critical transition, we know exactly where it is, and we know the
warning signal is present in the latent process. Can standard EWS methods find
it?**

If the answer were "easily", Phase 1 would be a formality. It is not.

---

## 1. The detector works on clean data

Gate passed (`sim/01`). On densely sampled, noise-free data the rolling
autocorrelation and variance both rise ahead of the transition, Kendall's tau
= 0.70-0.73 for every indicator, and the empirical rolling AC1 tracks the
analytic recovery rate `lambda(c)` derived from the model. On a non-tipping null
the same machinery gives tau = 0.01-0.03.

Euler-Maruyama step size was verified rather than assumed: dt = 0.01 reproduces
dt = 0.002 (mean transition day 304.6 vs 304.7) while dt = 0.02 visibly
disagrees (313.2).

So the pipeline is not broken. Everything that follows is about what reality
does to it.

---

## 2. The standard significance test is wrong by a factor of six

Rolling windows with a 2-day step and a 30-day width share ~93% of their data.
The Kendall tau computed on that indicator series is heavily autocorrelated, and
its nominal p-value assumes independence.

| | tau at the 5% level |
|---|---|
| nominal (analytic Kendall) | 0.12 - 0.14 |
| actual 5% point of the null | 0.36 - 0.67 |

**A nominally 5% test rejects 28-43% of the time under the null.** Every
threshold in this project is calibrated by simulation against a matched
no-change null, per design cell. Nominal p-values are banned (D9).

---

## 3. Measurement error creates a bias pointing the same way as the hypothesis

An AR(1) process observed with additive white noise is an ARMA(1,1); the naive
lag-1 autocorrelation is attenuated by the reliability ratio
`s2 / (s2 + sigma_m^2)`. Because `s2` grows as the system approaches the fold,
**the attenuation shrinks along the approach**, so the measured autocorrelation
rises partly for a reason that has nothing to do with slowing down.

Decomposed on a real trajectory (`sim/03`), observed AC1 rise +0.118:

- **+0.088** from genuine critical slowing down
- **+0.104** from the reliability ratio improving (0.25 -> 0.43)

**54% of the apparent "rising autocorrelation" is measurement error attenuating
less.** In a literature that tests the hypothesis by looking for a rise, this is
a bias in the direction of the conclusion.

---

## 4. The response scale can invent the signal, or reverse it

The largest finding in Phase 0, and it was found by chasing an anomaly: the
clinical null showed a *negative* variance trend (tau = -0.29) when the latent
variance in those same runs rises ~14%, exactly as theory predicts.

The variance of a rounded variable depends on where its mean sits between two
category boundaries — on a boundary you get near-Bernoulli variance across two
categories, mid-category you get almost none. **Approaching a transition the
latent mean slides, so it moves relative to the boundaries and drags rounding
variance with it**, by an amount that can dwarf and reverse the real signal.

`sim/07`, 200 replicates, identical dynamics throughout. Mean tau of the
variance indicator on the clinical null (correct answer: mildly positive):

| response model | tau on null | variance AUC |
|---|---|---|
| continuous | **+0.29** (correct) | 0.78 |
| Likert-7, euthymic mid-scale | +0.52 (inflated 1.8x) | 0.93 |
| Likert-7, euthymic near ceiling | **-0.33** (wrong sign) | 0.95 |
| Likert-11, euthymic near ceiling | **-0.63** (wrong sign) | 0.77 |

The autocorrelation indicators sit at AUC 0.72-0.76 across all four.

Confirmed inside the full design grid (`sim/04`, response model as an explicit
factor, mean AUC over 216 cells):

| indicator | continuous | Likert-7 | inflation |
|---|---|---|---|
| **variance** | 0.624 | 0.683 | **+0.058** |
| ac1_withinday | 0.595 | 0.595 | -0.001 |
| ac1_naive | 0.590 | 0.588 | -0.002 |
| ac1_ou | 0.598 | 0.592 | -0.006 |

Three consequences:

1. **Variance's apparent superiority is roughly tripled by the artefact.** Its
   honest lead over the best autocorrelation indicator is 0.624 vs 0.598 on
   continuous data, not 0.683 vs 0.592.
2. **A finer scale is not the fix.** Likert-11 was the worst configuration
   tested. What matters is where the mean sits and how far it travels relative
   to the boundaries, not how many boundaries there are.
3. **Rising variance on Likert ESM data is confounded with the mean's position
   on the response scale** — and a sliding mean is the defining feature of the
   situation the method is meant for.

---

## 5. Modelling the nuisance properly does not rescue it

`fit_ou_me` fits an OU state process with an explicit measurement-noise variance
by Kalman filter, identifying `lambda` separately from `sigma_m`. On parameter
recovery it does exactly what it should (`sim/05` Q1, true `lambda` = 0.6/hour):

| sigma_m | naive AC1 | `fit_ou` | `fit_ou_me` | sigma_m recovered |
|---|---|---|---|---|
| 0.1 | 0.34 | 0.78 | 0.62 | 0.08 |
| 0.4 | 0.10 | 3.36 | 1.13 | 0.32 |
| 0.8 | 0.03 | **12.45** | 3.70 | 0.68 |

`fit_ou` is catastrophic at realistic noise: it reads white measurement noise as
ultra-fast mean reversion and returns `lambda` 20x too large.

**And the corrected estimator still loses.** AUC against the clinical null:

| sigma_m | `ac1_ou` | `ac1_naive` | `ac1_ou_me` |
|---|---|---|---|
| 0.1 | 0.91 | 0.92 | **0.72** |
| 0.4 | 0.71 | 0.69 | **0.62** |
| 0.8 | 0.65 | 0.62 | **0.55** |

Worst at every noise level. Bias-variance: on rolling windows of ~450
observations a four-parameter Kalman fit is noisy enough that the indicator
series is mostly sampling error. Removing the bias cost more variance than the
bias was worth. This would presumably reverse with longer windows or estimation
pooled across people; neither is tested.

---

## 6. The power analysis: what is the minimum viable series?

216 design cells x 4 indicators x 3 scenarios x 400 replicates. Thresholds
calibrated per cell. "Viable" requires **both** 80% power and no more than 10%
false alarms on someone who deteriorates without becoming ill — a detector that
fires on everyone has power 1.

### The answer: there isn't one.

**No configuration in the grid clears both bars.** Not one.

Only two cells reach 80% power at all, and they buy it at a price no clinician
would pay:

| response | days | prompts/day | missing | sigma_m | indicator | power | **false alarms** |
|---|---|---|---|---|---|---|---|
| continuous | 180 | 6 | none | 0.1 | variance | 0.900 | **0.302** |
| continuous | 180 | 10 | none | 0.1 | variance | 0.877 | **0.260** |

Both require a continuous response, no missing data, and the lowest measurement
noise tested. Flagging 30% of people who were never going to become ill is not
a usable warning signal.

Conversely, every cell holding false alarms under 10% tops out at **0.68 power**
— and all of them are Likert-7, i.e. their low false-alarm rate is the
artefact of section 4 suppressing the null, not the detector's skill.

**171 of 216 design cells cannot reach even 30% power with any indicator.**

### What each design lever actually buys

| lever | | | | |
|---|---|---|---|---|
| days | 30: 0.075 | 60: 0.094 | 120: 0.177 | 180: 0.227 |
| prompts/day | 3: 0.101 | 6: 0.148 | 10: 0.183 | |
| measurement noise | 0.1: 0.223 | 0.4: 0.127 | 0.8: 0.087 | |
| missingness | none: 0.168 | MCAR 30%: 0.139 | MNAR 30%: 0.128 | |

(mean power over the rest of the grid)

**Measurement noise is the dominant lever** — a 2.6x swing, larger than
anything the study design controls. Reducing it means better instruments or
multi-item composites, not more beeps.

### The one clear design recommendation: buy duration, not density

Same prompt budget spent two ways (continuous response, no missingness, best
indicator):

| days | prompts/day | total prompts | power |
|---|---|---|---|
| **180** | **3** | **540** | **0.688** |
| 120 | 10 | 1200 | 0.482 |
| 120 | 6 | 720 | 0.515 |
| 180 | 6 | 1080 | 0.900 |

**180 days at 3 prompts/day beats 120 days at 10 prompts/day, using less than
half the prompts.** Critical slowing down is a slow phenomenon; you need to
observe the approach, and sampling the same short approach more finely does not
substitute for observing a longer one. This is the most actionable result in
Phase 0 and it is also the cheapest for participants.

---

## What this does and does not say

**It does not say the theory is wrong.** The planted signal is real and
detectable on clean data. Every failure below that is a measurement and
inference problem.

**It does say the standard analysis is far more fragile than the applied
literature suggests.** Three separate mechanisms — an anticonservative trend
test, a measurement-error bias pointing at the hypothesis, and a discretisation
artefact that can flip the sign of the primary indicator — all push in the
direction of finding critical slowing down whether or not it is there.

**These are best-case numbers.** The control parameter moves linearly and
monotonically toward the fold, which is the most favourable possible case for a
monotone trend statistic. Real deterioration is not monotone. The model is
one-dimensional, so it cannot produce the cross-correlation indicator used in
the published work at all. And the timescale anchor (D3) — one model time unit
= 5 hours, calibrated to an ESM autocorrelation of 0.4 — scales every statement
about days and is not swept.

**Implication for Phase 1.** The Kossakowski et al. data is a single person with
a single transition, ~10 prompts/day over 239 days, on Likert items, with
missingness. That sits close to the best cells in this grid — but those cells
had 26-30% false alarm rates, and the Likert cells' apparent specificity is
artefact. A successful reproduction there should be treated as consistent with
the hypothesis, not as evidence for it, and the specification curve in Phase 2
is where the actual test lives.

---

## Reproducing

```
Rscript sim/01_ideal.R              #  ~20 s   gate
Rscript sim/02_null_calibration.R   #  ~2 min  thresholds
Rscript sim/03_degrade.R            #  ~1 min  degradation + attenuation
Rscript sim/04_power.R              # ~11 min  the grid (cached; REFRESH=TRUE to redo)
Rscript sim/05_measurement_error.R  # ~28 min  state-space comparison
Rscript sim/07_likert_artifact.R    #  ~2 min  response-scale artefact
Rscript sim/06_summary.R            #   ~1 s   this document's numbers
```
