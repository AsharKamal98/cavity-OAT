#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Harmonic Analysis]
]

= Purpose

This file defines the solver-independent harmonic diagnostic implemented in
`post_analysis/harmonic_analysis.py` and its output containers in
`parser/harmonic_analysis.py`. Use it to analyze one or more clean real series,
including bright or dark mode components, on a shared time grid.

= Definitions

For one input series $y(t)$, fit an offset and integer harmonics of the
fundamental frequency $f_0$:

$
y(t) = c + sum_(n=1)^(N_h) A_n cos(2 pi n f_0 (t-t_0) + chi_n).
$

Here $t_0$ is the first supplied time, $A_n$ is the amplitude, and $chi_n$ is
stored as `phase_offsets`. Frequencies are stored in cycles per input time unit,
$f_n=n f_0$.

== Metrics of Interest

=== RMS Oscillation Amplitude

Define the RMS oscillation amplitude as

$
A_"RMS" = sqrt(frac(1, 2) sum_(n=1)^(N_h) A_n^2).
$

This measures the total strength of the fitted oscillatory part while excluding
the constant offset $c$. It combines harmonic amplitudes in quadrature, so it
does not depend on their fitted phases. A constant or numerically constant
series has $A_"RMS"=0$.

=== Total Harmonic Distortion

The total harmonic distortion (THD) is

$
"THD" = frac(sqrt(sum_(n=2)^(N_h) A_n^2), A_1).
$

THD measures the combined higher-harmonic amplitude relative to the fundamental
amplitude. A pure fitted sinusoid has $"THD"=0$; increasing THD means the fitted
waveform is less sinusoidal. THD is `NaN` when no fundamental oscillation
exists. The metric refers only to harmonics retained by `max_harmonic` and
supported below the Nyquist frequency.
If only the fundamental can be fitted, report THD as `NaN`, since no
higher-harmonic content was measured.

=== Fit Quality

Store `fit_r_squared` as one minus the squared fit residual divided by the
squared variation about the window mean. A value near one indicates a good
local periodic fit. This is a fit diagnostic, not a significance test or an
uncertainty estimate; correlated MCWF fluctuations can also look periodic.

= Method in Pseudo-code

```python
def extract_family_phase_series(t, series, phase_protocol, *,
                                family_phase_index):
    validate the zero-based family-phase index
    select the inclusive family-phase time interval
    return the selected time grid and selected curves

def compute_harmonic_analysis(t, series, *, max_harmonic=5):
    t, dt = validate finite, increasing, uniform sampling
    curves = as_series_tuple(series)

    analyses = []
    for curve in curves:
        f0 = estimate_fundamental_frequency(t, curve, dt, max_harmonic)
        if f0 == 0:
            append a non-oscillating result with empty harmonic arrays and
                THD = NaN
            continue
        if f0 * (t[-1] - t[0]) < 1:
            return None
        analyses.append(
            fit_harmonics(t, curve, f0, dt, max_harmonic)
        )

    return HarmonicAnalysisSeries(signals=tuple(analyses))

def estimate_fundamental_frequency(t, curve, dt, max_harmonic) -> float:
    centered = curve - mean(curve)
    fft = numpy.fft.rfft(centered)
    return 0 when the maximum centered amplitude is at most 1e-10 times
        max(1, maximum absolute curve value)
    frequency_guess = frequency of largest non-DC FFT magnitude
    refine frequency_guess between neighboring FFT bins by minimizing the
        least-squares residual of the requested harmonic model
    return refined_frequency

def fit_harmonics(t, curve, f0, dt, max_harmonic) -> HarmonicSignalAnalysis:
    retain integer harmonics below the Nyquist frequency
    fit offset plus cosine/sine coefficients simultaneously
    convert coefficients to amplitudes and phase offsets
    compute RMS oscillation amplitude
    compute total harmonic distortion
    return HarmonicSignalAnalysis(...)
```

`compute_harmonic_analysis(...)` normalizes the single-curve or multi-curve
input using `common.utils.moments.as_series_tuple(...)`. The FFT supplies only
the initial frequency estimate; continuous-frequency least-squares refinement
and the final simultaneous harmonic fit allow a partial final oscillation.

`extract_family_phase_series(...)` performs the shared phase-window selection
before analysis. `family_phase_index` is zero-based, so phase 2 uses index 1.
It returns `(phase_t, phase_series)` and keeps the supplied curve order.

== Windowed Analysis

`compute_windowed_harmonic_analysis(t, series, *, window_duration,
step_duration, max_harmonic=5, min_fit_r_squared=None)` repeatedly calls the existing
`compute_harmonic_analysis(...)` to track the offset, fundamental amplitude,
RMS amplitude, frequency, and THD over time. Durations use the same units as `t`.

```python
def compute_windowed_harmonic_analysis(t, series, *, window_duration,
                                     step_duration, max_harmonic=5,
                                     min_fit_r_squared=None):
    validate the shared time grid, curves, and positive finite durations
    construct complete windows starting at t[0], advancing by step_duration
    initialize each metric as NaN for every window and input curve
    for window in windows:
        if fewer than four saved samples are selected:
            continue
        for curve in curves:
            result = compute_harmonic_analysis(window.t, curve[window],
                                               max_harmonic=max_harmonic)
            if result is not None:
                record fit quality
                record offset, fundamental amplitude, RMS amplitude, frequency, THD
                    when the optional fit-quality threshold is satisfied
    return WindowedHarmonicAnalysisSeries(window_centers, metrics)
```

The wrapper owns window selection, fit-quality filtering, and collection of
fit results. Fit curves independently
so an unresolved oscillation leaves a gap only for that curve and window.
Use inclusive endpoints, allowing numerical roundoff at grid-aligned bounds;
omit an incomplete trailing window. The reported time is the requested window's
midpoint. A step shorter than the window duration produces overlapping windows.
With `min_fit_r_squared=None`, retain the original unresolved-window behavior.
With a threshold in $[0,1]$, retain the window mean as a background estimate
when the periodic fit is unresolved or rejected, and leave its oscillation
metrics as `NaN`. Also skip an unresolved Nyquist-edge peak in this mode.
Other invalid inputs still raise. Constant signals retain their established
zero-amplitude convention.

Select phase 2 with `extract_family_phase_series(...)` before calling the
wrapper when analyzing dark-mode decay. Each window fits a constant offset
and amplitude locally; these are estimates averaged over the window, rather
than instantaneous measurements. The offset is the sinusoidal midpoint; higher
harmonics can make it differ from the midpoint between extrema. Overlapping
window estimates are correlated.

`notebooks/bright_dark_modes_mcwf.ipynb` applies this wrapper to the three
phase-2 dark components for MCWF and MFE. Its SLURM counterpart,
`slurm/analyses/bright_dark_modes_mcwf/post_analysis.ipynb`, uses the same
windowed calls, quality flags, plots, and separate early-time bootstrap on
saved raw trajectories. Use the loaded artifact's physical metadata and time
grid, and check its resolution against the reference notebook's saved-time
spacing. Keep window and step durations
explicit next to the calls. Its validated example uses $T_2=5000$, saved
spacing at most $2/Gamma$, windows of $200/Gamma$, and steps of $50/Gamma$.
Derive the saved-point count from total duration so extending the run does not
silently coarsen the grid. These settings resolve the example's fundamental
near $0.0338 Gamma$; reassess them when physical parameters change.
Keep candidate fits with `min_fit_r_squared=0.0`. Separately flag fits below
0.9 R-squared, with fewer than three cycles per window, or with fewer than
eight samples per period. Print candidate counts, missing fits, each failure
count, accepted counts, and the R-squared range and median for every component.
Use `compute_harmonic_quality(analysis, sample_spacing=..., window_duration=...,
labels=..., verbose=True)` to compute the `(windows, curves)` boolean masks
and print these diagnostics. It accepts the three display thresholds as
optional arguments and does not modify or refit the stored metrics. Call it
explicitly for MCWF and MFE instead of repeating component loops in notebooks.
Failure counts can overlap. These editable criteria judge a local periodic
approximation; they do not establish whether a physical oscillation exists.
Reuse `plot_harmonic_sweep(...)` with window-centre time on its x-axis,
passing the boolean display criteria through `quality_masks`. Show accepted
metrics as lines and finite flagged candidates as crosses, without connecting
across gaps. Missing fits remain `NaN`. Always plot the background. Plot MCWF
solid and MFE dashed, with matching component colors, and explain the crosses
in the figure title. The shared plotter remains compatible with parameter
sweeps without a quality mask. Pass transposed metric and mask arrays to plot
all three components in one call per method. Plot $f_0/Gamma$ against $Gamma t$ in this
notebook. The background is the candidate fit offset, or the window mean when
no fit is available; label that fallback explicitly. A flagged amplitude or
frequency can reflect noise or model mismatch and must not be interpreted as
a resolved oscillation.

== Early-Time Damping and Trajectory Bootstrap

`compute_damped_oscillation(t, series, *, frequency_bounds)` fits one real
series to an exponential sinusoid with a linear background:

$
y(t) = c + b u + exp(-gamma u)
  (a cos(2 pi f u) + d sin(2 pi f u)), quad u=t-t_0.
$

The initial amplitude is $sqrt(a^2+d^2)$; the plotted envelope is this
amplitude times $exp(-gamma u)$ around the linear background. Frequencies are
cycles per input time unit. Fit the background and sine/cosine coefficients
by linear least squares for each trial decay rate and frequency, then choose
the nonlinear parameters by least squares. Allow signed decay rates so the
fit can return zero damping or growth. Limit the envelope's change over the
interval to a factor $exp(10)$ in either direction and report a
`at_parameter_bound` flag when that limit or a frequency bound is reached.
Numerically constant or purely linear series return `None`. Require a uniform
grid with at least eight samples and explicit positive frequency bounds below
Nyquist.

`compute_bootstrap_damped_oscillations(...)` computes the central fits and
refits after resampling complete independent trajectories with replacement.
Its callback receives trajectory indices, including repetitions, and returns
the desired curves. Preserve all times and groups for each sampled trajectory;
do not resample time points or treat overlapping windows as independent data.
The short-time analysis callback reuses `compute_mean_generator_moments(...)`,
then `compute_j_moments(...)`, then `compute_j_modes_from_j_moments(...)`.
Thus the original averaging and observable-construction order is preserved.

Keep the notebook orchestration in
`post_analysis/short_time_oscillations.py` through
`compute_short_time_dark_mode_analysis(moments, *, duration,
reference_frequencies, num_resamples=200, confidence_level=0.95, seed=None,
verbose=False)`. This helper selects the start of phase 2, caps the requested
duration at its end, prepares the cached trajectories, builds the resampling
callback, calls the existing bootstrap fitter, and optionally prints its
summary. It returns `ShortTimeDarkModeAnalysis` from
`parser/short_time_oscillations.py`: `t`, the dark x/y/z `components`, nullable
component `sem`, and the existing damped `fits`. Pass these fields directly to
`plot_damped_oscillations(...)`; keep data selection, callbacks, and reporting
out of notebook cells. Do not change the windowed analysis or generic fitter
to implement this orchestration.

The notebook requests the first $200/Gamma$ of cached phase-2 data, 200
resamples, 95% percentile intervals, and a reproducible seed. It supplies its
existing MFE frequencies as plain reference values. The helper sets the band
to half to one-and-a-half times the median positive reference frequency,
capped below Nyquist; it does not run or modify windowed analysis. Analyze
the cached raw moments without rerunning either simulation.
Never rerun parameter initialization merely to import an analysis helper:
that cell replaces the result containers.

Return amplitude, decay-rate, and frequency intervals only when the central
fit avoids parameter bounds and at least 80% of resampled fits converge away
from those bounds. Always report successful resample counts, fit R-squared,
and bound flags. A positive decay interval supports damping *within this
chosen model and band*; it does not establish that a single exponential fully
describes the dynamics. Inspect the fitted trace and residual mismatch before
interpreting it. Pointwise SEM bands in `plot_damped_oscillations(...)` describe
uncertainty in the MCWF mean, not the fitted envelope's confidence interval.
Pass cached MFE times and components through the plotter's `reference_t` and
`reference_series` arguments. It selects the matching fit interval and draws
the dashed reference in every panel; no notebook loop is needed.

The damped fit returns `DampedOscillationAnalysis`, including the time grid,
fitted curve, background, envelope, fitted parameters, fit quality, bound flag,
and nullable `DampedOscillationUncertainty`. The bootstrap returns these fits
in `DampedOscillationSeries.signals`, preserving curve order and using `None`
where no oscillatory fit is available. Keep these outputs separate from the
simulation and J-mode containers.

= Output

```python
HarmonicAnalysisSeries(
    signals=(
        HarmonicSignalAnalysis(
            offset,
            fundamental_frequency,
            observed_cycles,
            harmonic_indices,
            frequencies,
            amplitudes,
            rms_amplitude,
            phase_offsets,
            total_harmonic_distortion,
            fit_r_squared,
        ),
        ...,
    ),
)
```

`HarmonicAnalysisSeries` members are:

- `signals`: one `HarmonicSignalAnalysis` for each supplied series, in input
  order.

Each `HarmonicSignalAnalysis` contains:

- `offset`: the fitted constant background $c$.
- `fundamental_frequency`: the fitted base frequency $f_0$, in cycles per input
  time unit.
- `observed_cycles`: the approximate number of fundamental oscillations in the
  supplied interval, $f_0 (t_"last"-t_0)$.
- `harmonic_indices`: the retained harmonic numbers $n=1,2,dots$.
- `frequencies`: the corresponding frequencies $f_n=n f_0$.
- `amplitudes`: the fitted non-negative amplitudes $A_n$.
- `rms_amplitude`: the RMS oscillation amplitude $A_"RMS"$ defined above.
- `phase_offsets`: the fitted phases $chi_n$ relative to the first supplied
  time $t_0$.
- `total_harmonic_distortion`: the combined root-sum-square amplitude of all
  retained harmonics with $n >= 2$, divided by the fundamental amplitude. It
  is `NaN` when the series is constant because no fundamental exists.
- `fit_r_squared`: fraction of centered variation explained by the fit; one
  for the constant-signal convention.

The array fields `harmonic_indices`, `frequencies`, `amplitudes`, and
`phase_offsets` correspond element by element.

The result is plot-ready. Parameter-sweep arrays extracted from several results
may be passed to `common.plotting.harmonic_analysis.plot_harmonic_sweep(...)`.
Plotting functions must not repeat the FFT or harmonic fit.

For the sweep plot, use the stored `rms_amplitude` and `offset` directly.

The windowed wrapper returns a `WindowedHarmonicAnalysisSeries` from
`parser/harmonic_analysis.py`:

```python
WindowedHarmonicAnalysisSeries(
    window_centers,              # shape (windows,)
    offsets,                     # each metric: (windows, curves)
    fundamental_amplitudes,
    rms_amplitudes,
    fundamental_frequencies,
    total_harmonic_distortions,
    fit_r_squared,
)
```

Metric columns preserve input curve order, including a singleton column for
one curve. Unresolved fits or windows selecting fewer than four samples have
`NaN` metrics. Constant curves retain the fitted offset and have zero
fundamental amplitude, RMS amplitude, and frequency, with `NaN` THD. For example, plot
`result.offsets[:, 0]` and `result.fundamental_amplitudes[:, 0]` against
`result.window_centers` for a single supplied component.
The optional fit-quality threshold retains background means for rejected fits.
Pass the complete window grid and `NaN` metric gaps to the four-panel plotter.
`fit_r_squared` is nullable for compatibility with older constructed outputs;
the current wrapper always fills its array.

For a constant or numerically constant curve, the result stores
`fundamental_frequency=0`,
`observed_cycles=0`, `rms_amplitude=0`, empty harmonic arrays, and
`total_harmonic_distortion=NaN`. If any non-constant supplied curve contains
fewer than one fitted oscillation,
`compute_harmonic_analysis(...)` returns `None` for the complete multi-curve
analysis. A parameter-sweep caller should print or warn using its own parameter
context, such as the current `delta`, and skip that parameter value.

= Data Requirements

- One real one-dimensional series or a sequence of real series.
- One shared finite, strictly increasing, uniformly sampled time grid with at
  least four points.
- Every non-constant series must contain at least one observed oscillation;
  otherwise the result is `None`. Constant series use the non-oscillating
  output convention described above.

= Invariants and Edge Cases

- Harmonic analysis must remain independent of the simulation backend and may
  consume any compatible real series.
- The largest non-DC FFT magnitude is assumed to identify the fundamental
  neighborhood. This simple implementation therefore assumes the fundamental
  is the dominant Fourier component.
- Test Nyquist equality only for an even sample count; the last FFT bin for
  an odd count is below Nyquist and is usable.
- Aliasing cannot be detected or repaired reliably from a sampled trace alone.
  Check saved-time resolution against an independently resolved frequency;
  interpolation does not restore omitted oscillations. A noise-only trace can
  produce nonzero fitted amplitude and an arbitrary frequency. Use trajectory
  resampling or a deterministic reference before interpreting these as physics.
- The input interval need not contain an integer number of oscillations.
- Treat a curve as numerically constant when its maximum centered amplitude is
  at most `1e-10 * max(1, max(abs(curve)))`. This prevents solver roundoff from
  being identified as an oscillation.
- Phase offsets are defined relative to the first supplied time $t_0$.
- `frequencies` use cycles per input time unit.
- `max_harmonic` must be positive. Harmonics at or above the Nyquist frequency,
  or unsupported by the sample count, must not be fitted.
- Reject non-finite, nonuniformly sampled, or shape-mismatched input with a
  clear `ValueError`.
- Complex transverse signals are not supported by this real-series helper.
- Window duration must span at least three saved-time spacings and must not
  exceed the supplied time span; step duration must be positive. Invalid input
  or errors raised by the existing fitter still raise rather than becoming
  missing values. No automatic period choice or decay model is introduced.
