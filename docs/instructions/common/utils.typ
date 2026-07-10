#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Common Utility Helpers]
]

= Purpose

This file lists the shared helper functions in `common/utils/phases.py`,
`common/utils/parameters.py`, and `common/utils/moments.py`. Use it when
adding, moving, or reusing generic
helpers that should be available outside a specific simulation backend such as
`solvers/mcwf`.

Detailed physics conventions may live in more specific instruction files, such
as `docs/instructions/model_parameters.typ` or
`docs/instructions/bloch_vector_averaging.typ`.

= Helper Summary

== `common/utils/phases.py`

This file contains backend-neutral helpers used across phase handling.

=== Phase Helpers

- `default_three_phase_protocol(T1, T2, T3, delta0, Omega0)` returns the
  standard three-phase protocol as a list of `Phase` objects.
- `phase_boundary_times(phases)` returns all cumulative phase-end times for a
  phase protocol.
- `phase_values_at_time(t, phases)` returns the phase-local `(Omega, delta)`
  values for a piecewise-constant protocol.

== `common/utils/parameters.py`

This file contains shared parameter-scale, cavity-mapping, and
parameter-validation helpers.

=== Parameter Helpers

- `omega_c(N_J, Gamma)` returns the critical drive $Omega_c$ for a fixed active
  atom number.
- `scaled_N_Gamma(factor, N, Gamma)` returns `factor * N * Gamma` for
  direct $N$-scaling model parameters such as $Omega_0$ and $delta_0$.
- `inverse_scaled_N_Gamma(factor, N, Gamma)` returns `factor / (N * Gamma)`
  for inverse rate scales such as MCWF timesteps.
- `mcwf_dt_from_scales(Omega0, delta0, N, Gamma, ...)` returns the minimum of
  the drive-based, detuning-based, and collective-decay-based MCWF timestep
  scales.
- `Omega_Gamma_from_cavity_parameters(...)` converts cavity parameters to the
  effective spin-model `(Omega, Gamma)` and checks the bad-cavity condition.
- `omega2_from_weighted_average(omega1, N1, N2)` returns the fixed
  inhomogeneous group-2 coupling $omega_2$ from the weighted-average
  convention.
- `validated_mcwf_dt(dt, N, Gamma, safety_factor=...)` enforces the older
  simple `1 / (N Gamma)` notebook timestep bound.
- `check_initial_sector_omega_ratio(sector_coeffs, Omega, Gamma, ...)` checks
  the initial support against the smallest-sector critical drive.

== `common/utils/moments.py`

This file contains shared moment/vector conversion helpers used by the new
pipeline.

=== Spin-Conversion Helpers

- `norm_spin_components_from_spin_components(x, y, z, tol=...)` returns
  Euclidean vector length and normalized vector components.
- `angles_from_norm_spin_components(sx, sy, sz, valid, tol=...)` returns
  `(theta, phi)` from normalized spin components.

= Invariants

- Helpers in `common/utils/phases.py` should be backend-neutral and should not import
  from `solvers.mcwf` or `solvers.qutip_fixed_nj`.
- Helpers in `common/utils/parameters.py` should stay limited to shared
  parameter, protocol, and validation logic.
- Helpers in `common/utils/moments.py` should stay limited to shared
  moment/vector post-processing.
- This instruction file should stay organized by utility file, with one short
  subsection per file under `common/`.
- If a helper becomes backend-specific, move it out of `common/utils/phases.py` rather
  than hiding backend assumptions in the common layer.
- When adding a new public helper to `common/utils/phases.py`,
  `common/utils/parameters.py`, or `common/utils/moments.py`, add a short entry
  to this file.
