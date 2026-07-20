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

The array fields `harmonic_indices`, `frequencies`, `amplitudes`, and
`phase_offsets` correspond element by element.

The result is plot-ready. Parameter-sweep arrays extracted from several results
may be passed to `common.plotting.harmonic_analysis.plot_harmonic_sweep(...)`.
Plotting functions must not repeat the FFT or harmonic fit.

For the sweep plot, use the stored `rms_amplitude` and `offset` directly.

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
