# Krylov and Matrix-Exponential-Action Propagation

A blanket Krylov replacement would probably slow the current $N\le100$ runs.
A hybrid backend becomes useful for larger two-group sectors, roughly once the
sector dimension approaches $900$–$1000$.

The current solver already has two propagation paths:

- Full-`dt` steps use a precomputed matrix
  $U=\exp(-iH_{\mathrm{eff}}dt)$, then apply `U @ psi`
  ([`sim.py`](../../solvers/mcwf/sim.py#L154)).
- Partial steps and jump-time bisections use `expm_multiply`, which computes
  $\exp(A)\psi$ without constructing the full propagator
  ([`sim.py`](../../solvers/mcwf/sim.py#L137)). This is already a
  matrix-exponential-action method. [SciPy documents `expm_multiply`
  here](https://docs.scipy.org/doc/scipy-1.15.1/reference/generated/scipy.sparse.linalg.expm_multiply.html).

Representative benchmarks on this machine:

| Two-group size | Sector dimension | `U @ psi` | `expm_multiply` | Better |
|---:|---:|---:|---:|---|
| $N=12$ | 16 | 0.002 ms | 0.133 ms | Precomputed $U$ |
| $N=30$ | 72 | 0.004 ms | 0.150 ms | Precomputed $U$ |
| $N=100$ | 676 | 0.258 ms | 0.327 ms | Precomputed $U$ |
| $N=120$ | 961 | 0.444 ms | 0.359 ms | Matrix action |
| $N=150$ | 1482 | 1.116 ms | 0.466 ms | Matrix action |

At $N=100$, about 93% of the recorded propagation calls use the precomputed
path. Replacing all of those with `expm_multiply` would therefore likely make
propagation roughly 20–25% slower.

The advantage reverses for larger sectors because the sparse generator's
exponential becomes dense:

- Dimension 676: about 9 MB per propagator per phase.
- Dimension 961: about 18.5 MB.
- Dimension 1482: about 44 MB.

Those matrices are stored for every phase and sector
([`sim.py`](../../solvers/mcwf/sim.py#L254)), so ramps, additional sectors, and
multiprocessing can substantially increase memory use. This will matter even
more when the temporary one-pair sector restriction is removed
([`state_helpers.py`](../../solvers/mcwf/state_helpers.py#L254)).

Recommendation:

- Keep precomputed propagators for small and medium sectors.
- Use `expm_multiply` for large sectors, initially around dimension 900, with
  the threshold confirmed by whole-trajectory benchmarks.
- Do not initially use SciPy's explicit restarted Krylov routine. In spot tests
  it was slower than both existing methods. SciPy also notes that its
  convergence depends strongly on the spectrum. [SciPy Krylov
  documentation](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.sparse.linalg.funm_multiply_krylov.html).
- Implement this eventually as
  `propagator_backend="precomputed" | "action" | "auto"`, selected separately
  for each phase-sector pair.

So: useful for scaling beyond the present $N=100$ two-group calculation, but
not as a universal replacement for the existing precomputed propagators.
