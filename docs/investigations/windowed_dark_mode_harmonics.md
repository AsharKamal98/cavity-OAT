# Windowed dark-mode harmonic analysis

Investigation of `notebooks/bright_dark_modes_mcwf.ipynb` and
`post_analysis/harmonic_analysis.py`.

The original plot mixes genuine quantum relaxation with undersampling and
finite-ensemble fluctuations. Its nearly constant background is reasonable
after relaxation, but its wandering MCWF frequencies and small fitted
amplitudes do not establish persistent coherent oscillations.

## Parameters and reference calculation

The audit used the notebook's physical parameters: `Ni=(5,5)`, couplings
`(0.7,1.3)`, `Gamma=1`, `Omega0=1`, `delta0=0.2`, `dN=0`, and phase durations
`(10,5000,10)`. The MCWF timestep was `0.01`, with 64 trajectories, seed 1234,
8 processes, and saved spacing 2. This run took about 159 seconds. The existing
10,000-duration images were inspected, but their raw trajectory arrays were
not available; the new stochastic audit is therefore a separate 64-trajectory
run, not a reanalysis of the original 500 trajectories.

For an independent, noise-free quantum reference, evolve the density matrix
in each retained sector `(k,5-k)`, with probability
`comb(5,k)*comb(5,5-k)/comb(10,5)`, starting with no excited atoms. Use

\[
H=\Omega(0.7J_{1x}+1.3J_{2x})-\delta(N_{e1}+N_{e2}),\qquad
L=0.7J_{1-}+1.3J_{2-},
\]

and the Lindblad dissipator `Gamma * D[L]`. Apply phase 1 for 10 time units,
then evolve through phase 2. The six sector Hilbert-space dimensions sum to
56; the diagonal density-matrix blocks contain only 560 entries. Inter-sector
coherences are unnecessary for these J observables because their operators
preserve each group sector. Form the dark vector from the resulting ensemble
group moments using the notebook's active-population normalization.

The unshifted master-equation generator and the shifted MCWF generator agree
to `1.8e-15` in their largest matrix-element difference. Diagonalized
propagation agrees with direct matrix exponentials at selected times through
1000 to within `5.0e-12` in density-matrix entries. Thus the reference does not
depend on stochastic convergence or the MCWF integration timestep.

The physical model is documented in `docs/theory/main.tex`,
`docs/theory/appendix_cavity_model.tex`, and
`docs/theory/notes_inhomogeneous_coupling.tex`.

## What explains the observed curves?

### Frequency: the old grid aliases the oscillations

The dense-grid MFE fundamental is approximately `0.03380437 * Gamma`, with a
period of `29.58 / Gamma`. With 500 saved points over a total duration of
10,020, the saved spacing is 20.08016 and the Nyquist frequency is 0.0249002.
The true fundamental is above Nyquist. Its alias is

\[
f_{\rm alias}=\left|\frac{1}{\Delta t}-f_0\right|
\simeq 0.0159960\,\Gamma,
\]

which explains the roughly 0.016 frequency in the old plot. Increasing the
duration to 100,000 while retaining 500 points makes the aliasing worse.
Interpolation of old saved data cannot undo it. The new spacing 2 gives about
14.8 saved points per true fundamental period and resolves all five requested
harmonics below Nyquist for these parameters.

No noisy MCWF frequency should be interpreted as physical drift merely
because an FFT peak exists. A stationary fluctuating trace also has a largest
Fourier component, and a fit can jump between components or their aliases.

### Offset: stationary backgrounds and a compressed vertical scale

Every window independently fits its own offset; the code does not reuse or
freeze one coefficient across windows. A synthetic linearly drifting midpoint
is recovered with a maximum error of `7.4e-8` in the calibration below.

The exact quantum solution approaches a background near
`(0.04174, 0.05210, -0.26163)` by the end of phase 2. Small remaining changes
are real relaxation. Plotting x/y backgrounds around 0.04--0.05 together with
z around -0.26 visually compresses those changes. Separate component axes
make them easier to inspect. Longer windows also average away faster changes;
the fitted offset is a local background estimate, not an instantaneous value.

Bootstrap resampling of complete trajectories, preserving their time
correlations, puts the exact late-time backgrounds inside the 64-trajectory
run's approximate 95% background intervals. These intervals include uncertainty
in the sampled sector populations; saved times are not independent samples.

### MCWF amplitude: real early decay followed by a noise floor

The exact quantum mean has damped oscillations. Examples of contributing
Liouvillian modes are:

| Sector | Frequency / Gamma | Decay time * Gamma | Initial dark-x amplitude contribution |
|---|---:|---:|---:|
| (3,2) | 0.038492 | 31.4 | 0.01310 |
| (2,3) | 0.034161 | 109.2 | 0.003066 |
| (1,4) | 0.032926 | 974.2 | 0.0002479 |

The slowly surviving coherent contribution becomes much smaller than the
finite-trajectory fluctuations. Different decay rates also mean that a single
sinusoid with a constant amplitude is only a local approximation to the early
quantum transient.

As a noise control, split the 64 trajectories into two independent groups of
32, compute their ensemble dark vectors separately, and take half their
difference. This largely cancels the deterministic mean and has approximately
the same sampling-noise scale as the full 64-trajectory mean. For 500-duration
windows, the unfiltered fitted RMS amplitudes are:

| Component | Median RMS fitted to MCWF mean | Median RMS fitted to noise control |
|---|---:|---:|
| x | 0.00675 | 0.00744 |
| y | 0.00639 | 0.00667 |
| z | 0.00722 | 0.00622 |

Thus a nonzero, slightly varying fitted amplitude at late times can be a
positive fitting bias from sampling noise. It is not evidence that the ensemble
mean maintains an oscillation of that amplitude. More trajectories reduce
this floor; a longer trace alone does not remove it.
Here "noise" means sampling uncertainty in the estimated ensemble mean.
The underlying conditional trajectories have physical stochastic dynamics
and can retain oscillatory correlations. This control does not assert that
their fluctuations are white noise or that all spectral structure is unphysical.

### MFE amplitude: a small decrease survives convergence checks

For 1000-duration windows, the MFE dark-x RMS amplitude decreases from
0.007064289 to 0.007044688 between windows centered near 750 and 4500 after
phase 2 begins: approximately 0.277%. RK45 with `rtol=1e-11, atol=1e-13` and
DOP853 with `rtol=1e-12, atol=1e-14` agree in these amplitudes to about
`2e-12`. The fitted frequency remains near 0.03380437 in both solutions.
This small amplitude decrease is consistent with the MFE dynamics rather than
time-integration error. It should not be conflated with the much stronger
quantum damping.

MFE matches the quantum initial mean active populations, but does not retain
the correlations introduced by conditioning on total `N_J=5`. Agreement of
the quantum and MFE traces is therefore not required for this small system.

## Calibration with known signals

On a grid with spacing 2, using 500-duration windows and step 125:

| Signal | Observed accuracy |
|---|---|
| Constant sinusoid, f=0.0338 | Frequency error below 3e-10; offset error below 5e-11 |
| Exponentially decreasing amplitude, decay time 2000 | Decay recovered; local offset bias below 1.9e-5 |
| Midpoint increasing as 0.04 + 1e-6 t | Offset error below 7.4e-8 |
| Frequency increasing from 0.03025 to 0.03475 at window centers | Frequency error below 3e-6 |
| Noise with a fixed mean and no sinusoid | Still produces fitted amplitudes; fit R-squared only 0.03--0.13 |

When amplitude, midpoint, and frequency all vary simultaneously, a
200-duration constant-parameter fit has frequency bias up to about `2.1e-5`
and offset bias up to `2.1e-5` in the calibration. Windowed estimates are
averages over the window and remain model-dependent. Overlapping windows
are correlated, not independent measurements.

## Initial implementation changes

- The notebook now uses `T2=5000`, saved spacing at most `2/Gamma`, and
  200-duration windows with 50-duration steps. The point count follows the
  total duration instead of remaining fixed as the run is extended. Other
  physical parameters, 500 default trajectories, and 8 processes are retained.
- Both harmonic output containers expose fit R-squared. The window helper
  accepts an optional quality threshold. Rejected periodic fits leave
  frequency/amplitude/THD as NaN while retaining the window mean as a labeled
  background estimate. Without a threshold, the old unresolved-fit behavior
  is retained.
- The notebook requires R-squared at least 0.9, three cycles per window, and
  eight samples per fitted period before plotting oscillation metrics. No
  MCWF window in the 64-trajectory audit passes all criteria; all 97 MFE
  windows pass for each component. These criteria are conservative display
  filters, not confidence intervals or a general proof of an oscillation.
- The plotter preserves NaN gaps instead of connecting surviving windows.
  Backgrounds remain visible even where periodic metrics are unresolved.
- The last FFT bin of an odd-length window is correctly treated as below
  Nyquist. A true Nyquist-edge peak remains unresolved. THD is NaN when only
  the fundamental is supported; a previous zero there did not establish a
  sinusoidal physical signal.
- The notebook rejects analysis of old coarsely sampled trajectories and
  clears old outputs so stale aliased plots are not presented as new results.

The existing harmonic-analysis and plotting tests pass: 8 tests and 4 subtests.
The updated notebook's post-processing and plotting cells were also executed
against the actual 64-trajectory first moments, including its SEM pipeline.
No new pytest tests or benchmarks were added.

## Reanalysis of the cached 5,000-trajectory run

The user subsequently ran 5,000 trajectories with phase durations
`(10,1000,10)`, retaining the same physical parameters and saved spacing 2.
The revised analysis was executed against that live notebook cache. Neither
simulation nor parameter initialization was rerun, and the kernel was not
restarted. The original simulation containers and trajectory array objects
were preserved.

For the 17 complete 200-duration windows, x and y each produced 17 candidate
fits. Every candidate failed R-squared 0.9; cycle count and saved resolution
passed. R-squared min/median/max were `(0.267,0.431,0.744)` for x and
`(0.326,0.460,0.764)` for y. Dark z produced six candidate fits, all below
R-squared 0.9; two also had fewer than three cycles. Its other eleven windows
had no available oscillatory fit. All 17 MFE windows passed for each component.

The notebook now retains finite candidate metrics and shows flagged values
as crosses, accompanied by the rejection counts and R-squared summaries.
Thus an all-rejected set remains inspectable without presenting its metrics
as resolved physical oscillations. Missing candidates still leave gaps.

The first `200/Gamma` of phase 2 were also fitted to a single exponential
sinusoid with a linear background. The frequency band was
`(0.0169022,0.0507066)*Gamma`, based on the MFE frequency. A 200-resample
whole-trajectory bootstrap repeated the usual mean-generator-to-J-to-dark-mode
pipeline, with seed 1234:

| Component | Decay rate / Gamma | Frequency / Gamma | Fit R-squared | Usable resamples | 95% decay-rate interval / Gamma |
|---|---:|---:|---:|---:|---|
| x | 0.04655 | 0.03848 | 0.884 | 128/200 | Not reported |
| y | 0.02891 | 0.03769 | 0.853 | 200/200 | [0.01456, 0.03919] |
| z | 0.00296 | 0.03856 | 0.875 | 136/200 | Not reported |

Fits at a parameter bound are excluded from the intervals; at least 80% of
resamples must remain usable. Consequently, x and z do not get intervals.
The y interval supports positive early-time damping *within the chosen
single-exponential model and frequency band*. Its imperfect R-squared and
the visible residual structure warn against treating the fitted rate as a
unique physical decay constant. In particular, the exact quantum reference
above contains several decaying modes, and a changing background can dominate
z's total fit quality.

The new fit plot overlays the actual cached means, pointwise SEM bands,
fitted curves, linear backgrounds, and envelopes. These SEM bands are not
confidence bands for the fitted envelope. Known synthetic signals recover
positive decay, zero decay, and growth; a purely linear series returns no
oscillatory fit. The existing eight tests and four subtests still pass.

## Local artifacts

Figures are in `output/bright_dark_modes_mcwf/harmonic_audit/`:

- `quantum_reference.png`: early and late dense-grid traces against the exact
  quantum solution and MFE.
- `decay_and_noise_floor.png`: exact relaxation and the noise-control fit.
- `background_detail.png`: separate component axes for background changes.
- `validated_windowed_harmonics.png`: the updated notebook's windowed plot
  using the 64-trajectory audit, with unresolved oscillation metrics omitted.
- `dense_phase_2_traces.png`: the updated phase-2 trace plot with SEM bands.
- `flagged_windows_5000.png`: candidate fits and display flags from the
  cached 5,000-trajectory run.
- `early_decay_5000.png`: early-time fitted curves, backgrounds, envelopes,
  and MCWF SEM bands from that same cached run.

Raw audit arrays, spectral coefficients, calibration summaries, and the audit
scripts are kept locally under `/private/tmp/cavity_oat_harmonic_audit` and
`/private/tmp/cavity_oat_harmonic_*.py`. They are not production 500-trajectory
simulation artifacts and are not inserted into the notebook's saved outputs.
