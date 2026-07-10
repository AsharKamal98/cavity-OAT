# Post Analysis

This file does not need to follow the `agent-instruction-writer` skill format.
It should simply list functions in the `post_analysis/` package with short
summaries attached to them.

## `post_analysis/mfe_residuals.py`

- `compute_mfe_residuals(j_moments, *, parameters, tol=1e-12)` computes the
  two-group complex MFE residual series from already-averaged J moments and
  shared moment parameters. It returns an `MFEResidualSeries` when the input
  contains compatible two-group fields, otherwise `None`.

## `post_analysis/theory_benchmarks.py`

- `phase1_ss_angles_for_nj(Nj, Omega, Gamma)` returns the phase-1
  steady-state benchmark angles `(theta_ss, phi_ss)` for one active sector.
