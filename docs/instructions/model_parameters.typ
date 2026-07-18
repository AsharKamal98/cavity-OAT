#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Model Parameters: Implementation Instructions]]

= Purpose

This file describes how simulation parameters ($Omega, Gamma, delta$ etc.) should be constructed. It should be used when editing parameter helpers.

= Model Parameters Definitions

== N-Scaled Effective Parameters

When the effective parameters should scale with system-size $N$, the notebook should define dimensionless factors and convert them with
shared helpers from `common/utils/parameters.py`:

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

== TWA-Matched Parameters

For comparisons with the TWA cavity-plus-spin code, use these helpers from
`common/utils/parameters.py`:

```
Gamma_from_twa(N) -> Gamma
Omega_from_twa(N) -> Omega_eff
```

They fix $g_1=1$, $g_2=g_1/2$, $g_"eff"=(g_1+g_2)/2$, and
$alpha_"in"=1$, giving

$
Gamma = frac(4, 15 sqrt(N / 2)),
quad
Omega_"eff" = frac(4 sqrt(g_"eff") alpha_"in", sqrt(15 sqrt(N / 2))).
$

== Cavity-Derived Parameters

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

For theory, check against `docs/theory/appendix_cavity_model.tex`

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

- `scaled_N_Gamma` in `common/utils/parameters.py` returns a direct
  $N Gamma$-scaled model parameter from a dimensionless factor.
- `Gamma_from_twa(N)` and `Omega_from_twa(N)` return the fixed parameter
  convention used to match the TWA cavity-plus-spin simulation.
- `Omega_Gamma_from_cavity_parameters` in
  `common/utils/parameters.py` returns $(Omega,Gamma)$.
- `omega_c` in `common/utils/parameters.py` returns $Omega_c$.
- `phase1_ss_angles_for_nj` in `post_analysis/theory_benchmarks.py` returns
  $(theta_("ss"),phi_("ss"))$.
- `check_initial_sector_omega_ratio` in `common/utils/parameters.py`
  returns the drive-ratio validation dictionary.
- `mcwf_dt_from_scales` in `common/utils/parameters.py` returns the
  current MCWF timestep from the drive, detuning, and collective decay scales;
  see `docs/instructions/solvers/mcwf/method_parameters.typ`.
- `omega2_from_weighted_average` in
  `common/utils/parameters.py` returns $omega_2$.

= Invariants

- Parameter helpers should be reused instead of duplicating formulas in
  notebook cells or plotting utilities.
- `Gamma` should remain the canonical code name for $Gamma$.
- Solver parameter classes should receive physical model inputs through shared
  `SimulationMetadata`; method-specific controls such as MCWF `dt` remain
  explicit solver inputs.
- Inhomogeneous $omega_2$ should be fixed from physical group sizes, not
  from each instantaneous sector.
- If an existing helper almost fits a new parameter convention but would
  need changed semantics, pause and ask before repurposing it.
