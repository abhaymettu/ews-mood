# Decisions log

Every entry is a choice that a competent reviewer could have made differently.
Where a choice is swept rather than fixed, that is noted. Phase 2 turns most of
these into axes of a specification curve; Phase 0 fixes them and flags them.

Format: **ID — the choice**, then what else was available, then why this one, then
what it would take to change my mind.

---

## Phase 0 (simulation)

### D1 — Cubic double-well as the mood model
`dx = (-x^3 + x + c) dt + sigma dW`.

Alternatives: a two-dimensional system (positive and negative affect as coupled
variables, which is what van de Leemput et al. actually simulate); a
higher-order potential with more than two states; a regime-switching model with
no continuous state at all; an empirically-fitted model.

Why: it is the normal form of the fold bifurcation, so any smooth one-dimensional
system approaching a fold looks like this near the transition. It is the minimum
structure that produces the phenomenon we are testing detectors for.

What would change my mind: nothing about Phase 0's conclusions, which are about
detector behaviour near a fold and are generic in that neighbourhood. But
**this model is a caricature and cannot be used to argue that real mood is
bistable.** It assumes the conclusion of the substantive debate. Phase 0 answers
"if mood were like this, could we detect it", not "is mood like this".

A two-variable version is the obvious extension, because cross-correlation
between emotions is a third EWS indicator in the published work and a
one-dimensional model cannot produce it. Not built. This is the largest gap in
Phase 0 relative to the papers being reproduced.

### D2 — Linear ramp in the control parameter
Alternatives: a random walk in `c`; a ramp with a plateau; step changes at life
events.

Why: the CSD derivation assumes `c` moves slowly relative to the system's
relaxation time, and a linear ramp is the cleanest realisation. A random-walk
`c` would be more realistic and would substantially lower detectability, because
approach to the fold would no longer be monotone.

Consequence to keep in mind: **the linear ramp is the most favourable possible
case for a monotone trend statistic like Kendall's tau.** Real power is lower
than anything reported here.

### D3 — Time calibration: 1 model time unit = 5 hours
Anchored so that far from the bifurcation, the lag-1 autocorrelation at
90-minute spacing is ~0.4, matching typical ESM "emotional inertia" estimates.

Alternatives: anchoring on the observed duration of depressive episodes; on
recovery time after a daily stressor; not anchoring at all and reporting
everything in model units.

Why: every statement about "how many days of data you need" is meaningless
without this anchor, so it has to be made explicitly rather than left implicit.

**This is the single most consequential number in the study.** Reported
inertia estimates range from roughly 0.2 to 0.6 depending on item, sampling
interval and population. At 0.2 the system is fast relative to the sampling
grid and detection improves; at 0.6 it is slow and detection degrades. Not
swept in Phase 0. Should be.

### D4 — Detrending: Gaussian kernel, bandwidth scaled to series length
Bandwidth = max(7, days/8), so 22.5 days for a 180-day series and 7 days for a
30-day one.

Alternatives: no detrending; linear; fixed bandwidth; first differencing.

Why detrend at all: approaching the fold the equilibrium slides from x* ~ 1.15
toward x* ~ 0.58. Rolling variance on undetrended data picks up that slide and
reports "rising variance" with a completely static noise process. Without
detrending the test is a test of the trend in the mean, which is not CSD and
which any depression questionnaire would have caught anyway.

Why scaled rather than fixed: a fixed bandwidth confounds series length with
detrending strength. A 14-day bandwidth removes most of the signal from a
30-day series and almost none from a 180-day one, so "longer series detect
better" would partly be "longer series were detrended more gently".

Why this is contestable: detrending too hard removes the low-frequency variance
that rising variance consists of. There is no principled bandwidth. The
scaling constant (days/8) is arbitrary. **Phase 2 sweeps this and it is a prime
candidate for the finding being fragile.**

### D5 — Process noise and measurement noise are separate parameters
`sigma_p` acts inside the SDE; `sigma_m` is added to the response, is white, and
does not grow near the fold.

Why it matters more than it looks: an AR(1) signal observed with additive white
noise is an ARMA(1,1), and the naive lag-1 autocorrelation is attenuated by the
reliability ratio `s2/(s2 + sigma_m^2)`. Because `s2` grows as the system nears
the fold, that ratio *improves* along the approach, so the naive AC1 rises
partly for a reason that has nothing to do with slowing down.
`sim/03_degrade.R` measures the split: at a plausible noise level, **59% of the
apparent rise in autocorrelation is the attenuation shrinking rather than the
recovery rate slowing.**

This is a bias in the same direction as the hypothesis, in a literature that
tests the hypothesis by looking for a rise. It is the strongest argument in
Phase 0 for a state-space estimator (`ind_ou_me_ac1`) over a moment-based one.

Contestable: the size of `sigma_m` relative to `sigma_p`. It is swept (0.1,
0.4, 0.8 scale points), which is the right treatment, but the *plausible* value
for real ESM data is genuinely unknown and is not identified by a single
series without a measurement model.

### D6 — MCAR and MNAR calibrated to the same marginal missing rate
The MNAR mechanism is `P(miss) = plogis(a - 1.5 * x)`, with `a` solved so the
expected missing rate equals the MCAR rate.

Why: otherwise a comparison of mechanisms is confounded with sample size, and
"MNAR is worse" would partly mean "MNAR threw away more data".

Contestable: the slope 1.5, which sets how strongly missingness depends on
mood. Not swept in Phase 0.

### D7 — Two nulls, not one
- **static** (`c` constant): sets all thresholds. Isolates the purely
  statistical false-positive rate.
- **drift** (`c` descends, stops short of the fold): the reported false-alarm
  rate. A person who deteriorates and does not become ill.

Why two: the static null answers "is the test correctly sized", the drift null
answers "how often does a clinician get a false alarm". They differ by a large
factor, and reporting only the first would flatter the method substantially.

A constraint discovered rather than chosen: **the drift null cannot be placed
arbitrarily close to the fold.** At `c = -0.235` the Kramers escape rate is
0.15 per time unit, so a system parked there tips with near-certainty inside a
study-length horizon. "Sits just short of tipping for months" is not a state
that exists in this model. The null was moved to `c_trough = -0.035`
(`drift_gap = 0.35`), calibrated to 0 tips in 120 runs over 300 days. A harder
null would require conditioning on survival, which biases the false-alarm rate.

### D8 — Detectors hand-rolled, not taken from `earlywarnings`
`generic_ews()` requires an evenly-spaced `ts` object. ESM data is not evenly
spaced, and pretending otherwise is one of the artefacts being quantified. The
implementations here are base R and are checked against theory
(`sim/01_ideal.R` compares the rolling AC1 against the analytic `lambda(c)`).

Phase 1 will cross-check against `earlywarnings` on the evenly-spaced case,
where it is applicable. If our AC1 disagrees with theirs on data they can
handle, we are wrong.

### D9 — Empirical thresholds only; nominal p-values are banned
Rolling windows with a 2-day step and a 30-day width share ~93% of their data
with their neighbours. The Kendall tau computed on that series is heavily
autocorrelated and its nominal p-value is badly anticonservative.

Measured in `sim/02`: the nominal 5% critical value is tau ~ 0.12-0.14; the
actual 5% point of the null distribution is tau ~ 0.37-0.67 depending on window
width. **Using the nominal p-value gives a false-positive rate of 29-40%, not
5%.** Every threshold in this project is calibrated by simulation, per cell.

### D10 — Observation window ends AT the transition
For tipping runs, the analysed window is the D days immediately preceding the
tipping point.

Alternatives: a window ending some fixed lead time before the transition; a
prospective online alarm evaluated at every time point.

Why: it is the most favourable framing (all the data, right up to the event)
and it matches what the published retrospective analyses do, so Phase 1 is
comparable.

Its limitation is important: **it uses the transition time, which in a real
prospective setting you do not know.** It therefore answers "was the signal
there in hindsight", not "would you have been warned in time". The prospective
version, with lead time and alarms per person-month, is Phase 3.

### D11 — Euler-Maruyama at dt = 0.01 model time units
Justified by a convergence check in `sim/01_ideal.R` against dt = 0.002:
transition-time distribution and stationary SD agree to within Monte Carlo
error. Sampling is at ~3-minute resolution, ~30 substeps per ESM prompt.

### D12 — Likert response scale — CORRECTED, this was wrong
**Original decision:** fix K = 7 and don't sweep it, on the grounds that
`sim/03` showed rounding costs only tau 0.528 -> 0.468, small next to
measurement noise, and that ceiling effects were "mild by construction".

**That was wrong, and it is the most consequential error in Phase 0.**

Found by asking why the clinical null had a *negative* variance trend
(tau = -0.29) when the latent variance in those same runs rises ~14%, exactly
as critical slowing down predicts.

Mechanism: the variance of a rounded variable depends on where its mean sits
between two category boundaries. A mean sitting on a boundary produces
near-Bernoulli variance across the two categories; a mean mid-category produces
almost none. Approaching a transition the latent mean *slides*, so it moves
relative to the boundaries and drags rounding variance along with it — by an
amount that can dwarf, and reverse, the real change in process variance.

`sim/07`, 200 replicates, mean tau of the variance indicator on the clinical
null (correct answer: mildly positive):

| response model | tau on null | variance AUC |
|---|---|---|
| continuous | +0.29 (correct) | 0.78 |
| Likert-7, euthymic mid-scale | +0.52 (inflated) | 0.93 |
| Likert-7, euthymic near ceiling | **-0.33 (wrong sign)** | 0.95 |
| Likert-11, euthymic near ceiling | **-0.63 (wrong sign)** | 0.77 |

The autocorrelation indicators sit at AUC 0.72-0.76 across all four response
models. Variance ranges 0.77-0.95 on identical dynamics.

Three things follow:
1. **The `sim/04` conclusion that variance is the strongest indicator was
   contaminated.** The grid used one arbitrary scale placement, and it happened
   to be one where the artefact aligned with the signal. The grid has been
   re-run with the response model (continuous vs Likert-7) as an explicit
   factor so the artefact-free ranking is visible.
2. **More categories is not the fix.** Likert-11 was the worst configuration
   tested. What matters is where the mean sits and how far it travels relative
   to the boundaries, not how many boundaries there are.
3. **Rising variance on Likert ESM data is confounded with the mean's position
   on the response scale.** Since a sliding mean is the *defining* feature of an
   approach to a transition, this confound is present in exactly the situation
   the method is meant to be used in. Autocorrelation indicators are far less
   vulnerable and should be preferred on discretised data, despite being weaker
   on continuous data.

Still not swept: where the scale is centred, and whether composites of several
items (which raise effective resolution) rescue the variance indicator.

### D14 — Modelling measurement error correctly does not help
`ind_ou_me_ac1` fits an OU state process with an explicit measurement-noise
variance by Kalman filter, so lambda is identified separately from sigma_m. It
was built because D5 shows the naive AC1 is biased in the same direction as the
hypothesis being tested.

It works as designed on parameter recovery (`sim/05` Q1, stationary OU on the
real ESM sampling grid, true lambda = 0.6/hour):

| sigma_m | naive AC1 | lambda, `fit_ou` | lambda, `fit_ou_me` | sigma_m recovered |
|---|---|---|---|---|
| 0.1 | 0.34 | 0.78 | 0.62 | 0.08 |
| 0.4 | 0.10 | 3.36 | 1.13 | 0.32 |
| 0.8 | 0.03 | **12.45** | 3.70 | 0.68 |

`fit_ou` is catastrophic at realistic noise: it reads white measurement noise as
ultra-fast mean reversion and returns lambda 20x too large. The state-space
version recovers both parameters far better.

**And it still loses.** On discrimination (`sim/05` Q2), AUC against the
clinical null:

| sigma_m | `ac1_ou` | `ac1_naive` | `ac1_ou_me` |
|---|---|---|---|
| 0.1 | 0.91 | 0.92 | 0.72 |
| 0.4 | 0.71 | 0.69 | 0.62 |
| 0.8 | 0.65 | 0.62 | 0.55 |

Worst indicator at every noise level. The reason is bias-variance: on rolling
windows of ~450 observations, a four-parameter Kalman fit has enough sampling
variance that the resulting indicator series is mostly noise, and a Kendall tau
computed on it carries little information. Removing the bias cost more variance
than the bias was worth.

This is worth stating plainly because the instinct — "the naive estimator is
biased, so model the nuisance properly" — is correct in principle and wrong
here at ESM window sizes. It would presumably become right with longer windows
or pooled estimation across people. Not tested.

### D13 — A 120-day hold before the ramp, to avoid conditioning on late tippers
The control parameter is held at its starting value for 120 days before the ramp
begins, so the median transition falls near day 314 rather than day 190.

Why, and this one was a bug before it was a decision: the analysis window is the
D days ending at the transition. Without the hold, a 180-day window could only
be cut from paths that happened to tip after day 180 — about 70% of them — and
those are not a random 70%. They are the paths where noise did *not* tip the
system early, which means they travelled further along the ramp and slowed down
more before tipping. The 30-day cells had no such restriction. So "longer series
detect better" would have been part real and part an artefact of conditioning the
long cells on the easiest paths, and the comparison across series length is the
entire deliverable.

With the hold, the earliest transition across 200 paths is day 244, so every
length from 30 to 180 days has full history available on every path and no
conditioning happens at all. The ramp *rate* is unchanged (the hold is prepended,
not compressed into the same span), so the approach speed and hence the amount of
slowing down available to detect is exactly as before.

Cost: 40% longer simulations. Worth it.

---

## Open, not yet decided

- **Two-variable model** with cross-correlation between affects (D1). Needed to
  reproduce the third published indicator at all.
- **Sweeping the timescale anchor** (D3), which currently propagates unexamined
  into every "days needed" number.
- **Non-monotone control parameter** (D2), which is the realistic case and will
  lower power.
- **Multiple testing across indicators.** The grid reports the best indicator
  per configuration; if an analyst picks the best of four post hoc, the true
  false-positive rate is higher than the per-cell 5%. Not yet corrected.

---

## Errors made and corrected

Four design errors were made during this project and caught before they reached
a result. They are written up in full in `report/report.qmd` section
"Errors made, and what they would have cost", because three of the four would
have produced a wrong result that looked publishable.

| error | would have caused | caught by |
|---|---|---|
| clinical null placed just short of the fold (D7) | every false-alarm rate computed against a null that was actually a transition | assertion counting how many nulls tipped |
| no hold phase before the ramp (D13) | "longer series detect better" inflated by conditioning long windows on late-tipping paths | checking what fraction of paths each cell used |
| "ceiling effects are mild by construction" (D12) | recommending the variance indicator, then calling its artefact a real signal in the real data | chasing an implausible sign in the clinical null |
| degradation layers compared on one path | reporting that missing data improves detection | re-running across 103 paths before believing it |

Three of the four pointed toward finding critical slowing down. That direction
is not coincidence and is discussed in the report.
