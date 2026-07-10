#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1.0em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Method Parameters: Implementation Instructions]]

= Purpose

This file describes numerical method parameters for the simulation backends.
Use it when editing timestep choices, saved-time grids, solver tolerances, or
backend-specific run controls.

Physical and model parameters such as $Omega$, $Gamma$, $delta$, and coupling
choices live in `docs/instructions/model_parameters.typ`.

= Monte-Carlo Wave Function

== Timestep Validation for the N-Scaled Choice

- `num_snapshots`: Sets the number of saved times on the requested `t_eval` grid. It therefore controls how many time points are returned by the solver and later appear in post-processing outputs such as moment series.

- `dt`: for N-scaled step-size, use 
  ```
  mcwf_dt_from_scales(Omega0, delta0, N, Gamma,
                      drive_factor=0.01, decay_factor=0.1) -> dt
  ```

  The current implemented rule is

  $
  Delta t = min(
    frac(0.01, |Omega_0|),
    frac(0.01, |delta_0|),
    frac(0.1, N Gamma),
  ).
  $

  The first term resolves coherent drive evolution on the timescale
  $1 / |Omega_0|$. The second resolves coherent phase accumulation from the
  detuning term on the timescale $1 / |delta_0|$. The third resolves collective
  dissipative evolution on the timescale $1 / (N Gamma)$.

  In the current MCWF implementation, a jump is first detected by checking
  whether the non-Hermitian norm crosses the random threshold during one
  attempted step, and the jump time is then refined by a fixed ten-step
  bisection in `solvers/mcwf/sim.py`. That reduces the jump-time uncertainty
  from order `dt` to order `dt / 2^10`. Because the code localizes the crossing
  after detection, the coarser collective-decay prefactor `0.1 / (N Gamma)` is
  currently sufficient for the outer attempted step, while the bisection refines
  the actual jump time.

  This timestep rule should not be duplicated manually in notebook cells. If a
  debugging workflow deliberately bypasses or changes the rule, that choice
  should be local and explicit.
