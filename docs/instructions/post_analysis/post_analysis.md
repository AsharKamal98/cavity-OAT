# Post Analysis

This file does not need to follow the `agent-instruction-writer` skill format.
It should simply list functions in the `post_analysis/` package with short
summaries attached to them.

## `post_analysis/mfe_residuals.py`

- `compute_mfe_residuals(j_moments, *, parameters, tol=1e-12)` computes the
  two-group complex MFE residual series from already-averaged J moments and
  shared moment parameters. It returns an `MFEResidualSeries` when the input
  contains compatible two-group fields, otherwise `None`.

## `post_analysis/j_modes.py`

- `compute_j_modes(t, x_components, y_components, z_components, *,
  populations, omega_groups)` computes one- or two-group common, contrast,
  drive-bright, and drive-dark vectors from active-population-normalized group
  components and returns a `JModeSeries`.

## `post_analysis/harmonic_analysis.py`

- `compute_harmonic_analysis(t, series, *, max_harmonic=5)`
  accepts one real series or a sequence of real series on a shared uniform time
  grid. For each series it removes the mean for FFT frequency discovery,
  refines the fundamental frequency between FFT bins, fits an offset plus
  integer harmonics, and returns a `HarmonicAnalysisSeries`.
- The fitted output stores harmonic indices, frequencies, amplitudes, RMS
  oscillation amplitude, phase offsets relative to the first supplied time,
  and total harmonic distortion.
  Constant or numerically constant series use fundamental frequency zero,
  empty harmonic arrays, and `NaN` THD. A non-constant series requires at
  least one observed oscillation; otherwise the function returns `None`. The
  interval need not contain an integer number of oscillations.

## `post_analysis/theory_benchmarks.py`

- `phase1_ss_angles_for_nj(Nj, Omega, Gamma)` returns the phase-1
  steady-state benchmark angles `(theta_ss, phi_ss)` for one active sector.
