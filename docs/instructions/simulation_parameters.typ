#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Simulation Parameters: Implementation Instructions]]

= Scope

This file describes how simulation parameters should be constructed in the current implementation. It should be used when editing notebook-level run functions, parameter helpers, phase construction, cavity-parameter conversion, or validation logic.

For physics derivations, use the files in `docs/theory`. In particular, cavity-derived parameters should be checked against `docs/theory/appendix_cavity_model.tex`. This file only states how those quantities should be used by the code.

= Overall Data Flow

Parameter setup should happen before constructing the initial sector state or starting a long trajectory run:

```
N, Gamma, dt, timing inputs, drive/detuning inputs
    -> Omega0, delta0
    -> phases
    -> validation checks
    -> run_trajectory_ensemble(...)
```

Notebook functions should take important run parameters, such as `N`, `Gamma`, and `dt`, as explicit inputs instead of relying on hidden global variables.

If the effective model parameters are chosen directly, the run function should pass $Omega$, $Gamma$, and $delta$ explicitly through to phase construction and simulation helpers.

= Effective Spin Parameters

For direct effective-model scans that use the current $N$-scaling convention,
the notebook should define dimensionless factors and convert them with shared
helpers from `common/utils/parameters.py`:

$
Omega_0 = c_(Omega) thin N Gamma, delta_0 = c_(delta) thin N Gamma.
$

The shared scaling helper is:

```
scaled_N_Gamma(factor, N, Gamma) -> value
```

so notebook cells should use:

```
Omega_factor, delta_factor
    -> Omega0 = scaled_N_Gamma(Omega_factor, N, Gamma)
    -> delta0 = scaled_N_Gamma(delta_factor, N, Gamma)
```

The older benchmark helpers

```
delta0_from_N_Gamma(N, Gamma) -> delta0
Omega0_from_N_Gamma(N, Gamma) -> Omega0
```

may still be reused when a notebook intentionally wants those fixed benchmark
choices. If a notebook cell instead chooses $Omega_0$ and $delta_0$ from
explicit scan factors, that choice should stay explicit and should not be
silently replaced by the older benchmark helper.

= Cavity-Derived Parameters

When the effective parameters should be derived from cavity inputs, code should reuse:

```
Omega_Gamma_from_cavity_parameters(
    epsilon, g_c, kappa, N_J, delta=0.0,
    bad_cavity_factor=10.0, round_digits=6,
) -> (Omega, Gamma)
```

This helper currently computes

$
Omega = frac(4 epsilon g_c, kappa), Gamma = frac(4 g_c^2, kappa).
$

It also checks the bad-cavity condition used by the implementation:

$
kappa >= mono("bad_cavity_factor") max(sqrt(N_J) thin |g_c|, |delta|).
$

If the condition is not satisfied, the helper should print the relevant values and stop the notebook-style workflow by raising `SystemExit(1)`. If the condition is satisfied, it should print the rounded $Omega$, $Gamma$, and bad-cavity diagnostic values and return $(Omega,Gamma)$.

= Critical Drive and Steady-State Helpers

The critical drive for a homogeneous active sector should be computed with:

```
omega_c(N_J, Gamma) -> Omega_c
```

where

$
Omega_c(N_J) = frac(1, 2) N_J Gamma.
$

The phase-1 steady-state angle diagnostic should reuse:

```
phase1_ss_angles_for_nj(Nj, Omega, Gamma) -> (theta_ss, phi_ss)
```

For $|Omega/Omega_c|<= 1$, this helper returns

$
theta_("ss") = arccos(sqrt(1-(Omega/Omega_c)^2)), phi_("ss") = frac(pi, 2).
$

For an invalid ratio, it returns `nan` values rather than pretending the steady-state angle exists.

= Initial-Sector Drive Validation

Before long scans over sector windows, code should use:

```
check_initial_sector_omega_ratio(sector_coeffs, Omega, Gamma,
                                 ratio_limit=1.0) -> dict
```

The helper should find the smallest populated total active number in `sector_coeffs`. For homogeneous keys this is $N_J$. For inhomogeneous keys $(N_(J,1),N_(J,2))$, this is $N_(J,1)+N_(J,2)$. It then evaluates

$
frac(Omega, Omega_c(N_(J,min))).
$

The returned dictionary should include:

```
is_valid, min_nj, omega, omega_c, ratio, ratio_limit
```

Notebook workflows may stop before entering a scan if `is_valid` is false.

= MCWF Timestep Validation

For custom MCWF runs using the current notebook timestep rule, code should
reuse:

```
mcwf_dt_from_scales(Omega0, delta0, N, Gamma,
                    drive_factor=0.01, decay_factor=0.1) -> dt
```

The current implemented rule is

$
dt = min(
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

= Phase Protocol

The standard three-phase protocol should be built with:

```
default_three_phase_protocol(T1, T2, T3, delta0, Omega0) -> phases
```

It returns a list of `Phase` objects from `parser/common.py`. Each `Phase` has:

```
Phase(duration: float, omega: float, delta: float, label: str = "")
```

The standard protocol is:

```
phase1: duration=T1, omega=Omega0, delta=0
phase2: duration=T2, omega=Omega0, delta=delta0
phase3: duration=T3, omega=0,      delta=0
```

Phase-boundary plotting or diagnostics should reuse:

```
phase_boundary_times(phases) -> [t1, t2, ...]
```

instead of recomputing the cumulative times manually.

Phase-local parameter lookup should reuse:

```
phase_values_at_time(t, phases) -> (Omega_t, delta_t)
```

instead of duplicating phase-boundary comparisons.

= Inhomogeneous Coupling Parameters

For two-group inhomogeneous coupling, group 2's coupling should be chosen with:

```
omega2_from_weighted_average(omega1, N1, N2) -> omega2
```

The current convention fixes the physical atom-number weighted average coupling to one:

$
N_1 omega_1 + N_2 omega_2 = N_1 + N_2.
$

This gives one fixed $omega_2$ for the whole inhomogeneous run. It should not be recomputed separately for each group-resolved sector. If $N_2=0$, the helper returns $1.0$ so downstream metadata remains finite.

For the inhomogeneous Hamiltonian, jump operator, sector keys, and residual diagnostics, use the dedicated inhomogeneous instruction and theory files rather than re-deriving those rules here.

= Reusable Helper Summary

- `delta0_from_N_Gamma` in `common/utils/parameters.py` returns
  $delta_0$.
- `Omega0_from_N_Gamma` in `common/utils/parameters.py` returns
  $Omega_0$.
- `scaled_N_Gamma` in `common/utils/parameters.py` returns a direct
  $N Gamma$-scaled model parameter from a dimensionless factor.
- `Omega_Gamma_from_cavity_parameters` in
  `common/utils/parameters.py` returns $(Omega,Gamma)$.
- `omega_c` in `common/utils/parameters.py` returns $Omega_c$.
- `phase1_ss_angles_for_nj` in `post_analysis/theory_benchmarks.py` returns
  $(theta_("ss"),phi_("ss"))$.
- `check_initial_sector_omega_ratio` in `common/utils/parameters.py`
  returns the drive-ratio validation dictionary.
- `mcwf_dt_from_scales` in `common/utils/parameters.py` returns the
  current notebook MCWF timestep from the drive, detuning, and collective
  decay scales.
- `default_three_phase_protocol` in `common/utils/phases.py`
  returns the standard phase list.
- `phase_boundary_times` in `common/utils/phases.py` returns all
  cumulative phase-end times.
- `phase_values_at_time` in `common/utils/phases.py` returns phase-local
  $(Omega(t),delta(t))$.
- `omega2_from_weighted_average` in
  `common/utils/parameters.py` returns $omega_2$.

= Invariants

- Parameter helpers should be reused instead of duplicating formulas in
  notebook cells or plotting utilities.
- `Gamma` should remain the canonical code name for $Gamma$.
- Functions that run simulations should receive `N`, `Gamma`,
  and `dt` explicitly when those quantities affect the result.
- The standard phase protocol should be represented with `Phase`
  objects, not with unrelated dictionaries or ad-hoc tuples.
- Inhomogeneous $omega_2$ should be fixed from physical group sizes, not
  from each instantaneous sector.
- If an existing helper almost fits a new parameter convention but would
  need changed semantics, pause and ask before repurposing it.
