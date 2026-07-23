# Dyadic Bisection Propagators

Dyadic propagators are worthwhile for jump-heavy trajectories with small or
medium sectors, but should not be a global default. In the current jump-heavy
benchmark they reduce trajectory runtime by about 49%; at sector dimension 676,
the setup and memory costs outweigh the much smaller propagation benefit.

## Method

The current four-iteration bisection performs five `expm_multiply` calls per
jump: four midpoint trials and one final propagation to the selected upper
bound. Instead, precompute

\[
U_k=\exp[-iH_{\mathrm{eff}}dt/2^k],\qquad k=1,\ldots,4.
\]

For a full-`dt` attempted step, retain states at the lower and upper bisection
bounds. At iteration $k$, apply $U_k$ once to the lower-bound state. The
result becomes either the new lower state or the new upper state. The final
upper state is already available, reducing each eligible localization from five
exponential actions to four matrix-vector products. Non-`dt` attempted steps
must retain the current variable-step path.

## Results

Measurements are medians on local arm64 macOS, Python 3.12.3 and SciPy 1.17.1,
using three integration phases and four bisection iterations. Memory is the CSC
storage of propagator matrices only.

| System | Sectors | Max dimension | Normal precompute | Extra dyadic, direct | Extra dyadic, square from $U_4$ | Existing memory | Extra memory |
|---|---:|---:|---:|---:|---:|---:|---:|
| $N=16$ | 3 | 30 | 10.8 ms | 36.8 ms | 9.9 ms | 0.08 MiB | 0.32 MiB |
| $N=40$ | 5 | 144 | 79.1 ms | 277.4 ms | 113.5 ms | 2.95 MiB | 11.80 MiB |
| $N=100$ | 1 | 676 | 314.6 ms | 847.6 ms | 1585.9 ms | 17.66 MiB | 70.65 MiB |

Computing $U(dt/16)$ once and squaring upward is faster for the two smaller
cases. Direct exponentiation is faster at dimension 676. The maximum relative
matrix difference between direct construction and squaring was
$2.2\times10^{-15}$.

| System | Current bisection | Dyadic bisection | Speedup | Eligible localizations to repay fastest setup |
|---|---:|---:|---:|---:|
| $N=16$ | 2.179 ms | 0.023 ms | 93.4x | 5 |
| $N=40$ | 4.181 ms | 0.184 ms | 22.8x | 29 |
| $N=100$ | 1.297 ms | 0.885 ms | 1.47x | 2058 |

The $N=40$ case has five sector blocks while $N=100$ has one, so these rows
should be interpreted as complete populated-sector cases rather than pure
dimension scaling.

The fixed-seed jump-heavy benchmark gives the end-to-end trajectory result:

| Path | Runtime | Jumps | `expm_multiply` calls |
|---|---:|---:|---:|
| Current | 86.3 ms | 20 | 177 |
| Dyadic prototype | 43.7 ms | 20 | 87 |

Eighteen of the twenty jumps occurred during full-`dt` attempts and used the
dyadic path; two occurred during partial steps and used the existing fallback.
Jump times were identical, and the final-state L2 difference was
$2.9\times10^{-15}$.

## Recommendation

- Implement dyadic bisection as opt-in first, or automatically only for
  phase-sector collections with maximum dimension around 150 or below.
- For small sectors, construct $U(dt/16)$ and square upward; retain direct
  exponentiation as the candidate for larger sectors.
- Use the incremental lower/upper-state algorithm so no final propagation is
  required.
- Fall back to `expm_multiply` whenever the attempted step is not exactly
  `dt`.
- Do not enable this for large sectors without an ensemble-level benchmark.
  Four dyadic levels add four times the existing propagator memory, potentially
  per worker process.

No solver code was changed during this investigation.
