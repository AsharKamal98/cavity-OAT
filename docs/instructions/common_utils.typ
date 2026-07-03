#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Common Utility Helpers]
]

= Purpose

This file lists the shared helper functions in `common/utils.py`,
`common/utils_parameters.py`, and `common/utils_moments.py`. Use it when
adding, moving, or reusing generic
helpers that should be available outside a specific simulation backend such as
`quantum_trajectories`.

Detailed physics conventions may live in more specific instruction files, such
as `docs/instructions/simulation_parameters.typ` or
`docs/instructions/bloch_vector_averaging.typ`.

= Helper Summary

== `common/utils.py`

This file contains backend-neutral general helpers used across phase handling,
legacy angle conventions, and observable comparisons.

=== Phase Helpers

- `phase_change_times(phases)` returns the first two phase-boundary times for a
  phase protocol.
- `phase_values_at_time(t, phases)` returns the phase-local `(Omega, delta)`
  values for a piecewise-constant protocol.

=== Legacy Angle Helpers

- `phase1_ss_angles_for_nj(Nj, Omega, Gamma)` returns the phase-1 steady-state
  angles `(theta_ss, phi_ss)` for one active sector.
- `active_manifold_angles(Jx, Jy, Jz, N_e, tol=...)` converts active-manifold
  observables into angles, active population, and normalized components.

== `common/utils_parameters.py`

This file contains shared parameter-scale, phase-protocol, cavity-mapping, and
parameter-validation helpers.

=== Parameter and Protocol Helpers

- `default_three_phase_protocol(T1, T2, T3, delta0, Omega0)` returns the
  standard three-phase protocol as a list of `Phase` objects.
- `omega_c(N_J, Gamma)` returns the critical drive $Omega_c$ for a fixed active
  atom number.
- `delta0_from_N_Gamma(N, Gamma)` returns the default notebook detuning scale.
- `Omega0_from_N_Gamma(N, Gamma)` returns the default notebook drive scale.
- `Omega_Gamma_from_cavity_parameters(...)` converts cavity parameters to the
  effective spin-model `(Omega, Gamma)` and checks the bad-cavity condition.
- `omega2_from_weighted_average(omega1, N1, N2)` returns the fixed
  inhomogeneous group-2 coupling $omega_2$ from the weighted-average
  convention.
- `validated_mcwf_dt(dt, N, Gamma, safety_factor=...)` enforces the notebook
  MCWF timestep rule.
- `check_initial_sector_omega_ratio(sector_coeffs, Omega, Gamma, ...)` checks
  the initial support against the smallest-sector critical drive.

== `common/utils_moments.py`

This file contains shared moment/vector conversion helpers used by the new
pipeline.

=== Spin-Conversion Helpers

- `norm_spin_components_from_spin_components(x, y, z, tol=...)` returns
  Euclidean vector length and normalized vector components.
- `angles_from_norm_spin_components(sx, sy, sz, valid, tol=...)` returns
  `(theta, phi)` from normalized spin components.

= Invariants

- Helpers in `common/utils.py` should be backend-neutral and should not import
  from `quantum_trajectories` or `quantum_trajectories_qutip`.
- Helpers in `common/utils_parameters.py` should stay limited to shared
  parameter, protocol, and validation logic.
- Helpers in `common/utils_moments.py` should stay limited to shared
  moment/vector post-processing.
- This instruction file should stay organized by utility file, with one short
  subsection per file under `common/`.
- If a helper becomes backend-specific, move it out of `common/utils.py` rather
  than hiding backend assumptions in the common layer.
- When adding a new public helper to `common/utils.py`,
  `common/utils_parameters.py`, or `common/utils_moments.py`, add a short entry
  to this file.
