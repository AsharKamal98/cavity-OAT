#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[J-Vector Modes]
]

= Purpose

This file defines the one- and two-group common, contrast, drive-bright, and drive-dark
vectors and their available SEMs computed in `post_analysis/j_modes.py`.
Use these derived modes for post-simulation comparison and plotting without
modifying raw `moments.J` data.

Prefer `compute_j_modes_from_j_moments(...)` for shared J-series inputs.
The existing array-based `compute_j_modes(...)` accepts already-normalized
vectors and returns central values without sampling uncertainty.

= Definitions

The inputs are the active-population-normalized group vectors

$
bold(s)_g(t) = frac(2 bold(J)_g(t), N_(J,g)(t)).
$

For two groups with active populations $N_(J,1)$ and $N_(J,2)$, define

$
bold(s)_"common" = frac(N_(J,1) bold(s)_1 + N_(J,2) bold(s)_2, N_(J,1) + N_(J,2)),
quad
bold(s)_"contrast" = frac(N_(J,2) bold(s)_1 - N_(J,1) bold(s)_2, N_(J,1) + N_(J,2)).
$

For coupling-dependent modes, let $w_g = N_(J,g) omega_g$ and
$W = abs(w_1) + abs(w_2)$. Then

$
bold(s)_"bright" = frac(w_1 bold(s)_1 + w_2 bold(s)_2, W),
quad
bold(s)_"dark" = frac(w_2 bold(s)_1 - w_1 bold(s)_2, W).
$

The dark coefficients $(w_2, -w_1)$ are orthogonal to the bright coefficients
$(w_1, w_2)$ in two-group coefficient space. For every mode $m$, the common
`BlochVectorSeries` container derives

$
L_m(t) = norm(bold(s)_m(t)),
quad
L_("xy",m)(t) = sqrt(s_(m,x)(t)^2 + s_(m,y)(t)^2).
$

For one group, $bold(s)_"common"=bold(s)_1$,
$bold(s)_"bright"="sign"(omega_1) bold(s)_1$, and the contrast and dark
vectors are zero. Common and bright coincide for positive coupling.

== Sampling SEMs

The J-series converter takes raw `x/y/z` and `N_j` group fields and performs
the $2 bold(J)_g/N_(J,g)$ normalization internally. All populations must
exceed `tol`; normalized modes are undefined at zero active population.

When `JMomentSeries.uncertainty` is present, use its full `(time, 5G, 5G)`
sample covariance in group-major `(x, y, z, N_e, N_j)` order. For each mode,
construct the Cartesian Jacobian with respect to these source fields. Include
derivatives of both the normalized inputs and the population-dependent mode
weights. Couplings are fixed metadata. Retain population/component and
cross-group correlations in propagation.

Follow the generic covariance and SEM flow in
`docs/instructions/post_analysis/uncertainty.typ`. Cartesian mode propagation
is exact when populations are fixed; otherwise it is first-order. Pass the
propagated `(time, 3, 3)` sample covariance and original sample count to each
mode's `BlochVectorSeries.uncertainty`. That container owns component SEMs and
first-order length and xy-length SEMs, using the shared `vector_length_sem`
helper in `common/utils/moments.py`. Undefined norm derivatives give `NaN`
SEMs at zero length. For a single group the zero contrast/dark components have
zero SEMs, while their norm SEMs follow this same undefined-derivative policy.

Missing source covariance leaves every mode uncertainty and SEM unset,
including MFE output. Never use the input SEM arrays as substitutes for full
covariance. Each mode retains its own Cartesian covariance; a later quantity
combining modes must return to the original J covariance and combined
Jacobians to preserve between-mode correlations.

= Method in Pseudo-code

```python
def compute_j_modes(t, x_components, y_components, z_components,
                    *, populations, omega_groups) -> JModeSeries:
    require matching x/y/z component series for one or two groups
    common = population-weighted average of the group vectors
    contrast = population-weighted relative difference of the group vectors
    bright, dark = coupling-weighted orthogonal combinations
    return JModeSeries(t, common, contrast, bright, dark)

def compute_j_modes_from_j_moments(j_moments: JMomentSeries,
                                 *, omega_groups, tol=1e-12) -> JModeSeries:
    require one or two groups and finite positive active populations
    normalize raw Cartesian group vectors by their active populations
    modes = compute_j_modes(t, normalized components,
                            populations=N_j, omega_groups=omega_groups)
    if source sampling covariance is present:
        construct each mode Jacobian including population derivatives
        propagate the full source covariance with the shared generic helper
        reconstruct each BlochVectorSeries with its covariance and sample count
            # The vector container derives component and norm SEMs.
    return modes
```

`compute_j_modes_from_j_moments(...)` owns normalization and uncertainty
propagation, reusing `compute_j_modes(...)` for central values. Both return
all four modes as one `JModeSeries`. The shared vector container owns length
and SEM attachment; parser details live in `docs/instructions/parser.typ`.

= Output

```python
JModeSeries(
    t,
    common=BlochVectorSeries(x, y, z, length, xy_length),
    contrast=BlochVectorSeries(x, y, z, length, xy_length),
    bright=BlochVectorSeries(x, y, z, length, xy_length),
    dark=BlochVectorSeries(x, y, z, length, xy_length),
)
```

Store this result as `moments.J_modes`; keep `moments.J` as the raw shared
J-moment representation. Each nested vector also has nullable `uncertainty`,
`x_sem`, `y_sem`, `z_sem`, `length_sem`, and `xy_length_sem` fields.

= Data Requirements

- One shared time grid and matching normalized x/y/z component series for one
  or two groups. Each component input may be one series or a collection of
  group series.
- One active-population series and completed coupling per physical group.
- For the J-series converter, matching raw group x/y/z fields and optional
  full J sampling covariance plus sample count.

= SLURM Analysis

`slurm/analyses/bright_dark_modes_mcwf/submit.slurm` runs one MCWF ensemble
through the shared `slurm/run_solvers/run_mcwf.py` entry point. Follow the
`squeezing_nj_ratio-scan` artifact workflow, but use a single run with index
zero rather than a ratio sweep. Save trajectory-resolved moments with
`--save-ensemble` to the analysis-local `output` directory; the filename is
`<Ni><job-name>_0.pkl`. The input notebook must use the same prefix and job
name as the submission script. Keep the physical inputs, initial-sector
filters, trajectory count, and seed aligned with
`notebooks/bright_dark_modes_mcwf.ipynb`; omitted optional filters remain
unset. The cluster worker count follows `SLURM_CPUS_PER_TASK`. Derive the saved
point count from the total duration and the notebook's saved-time spacing,
rather than fixing it at 500. Request memory for the full trajectory-resolved
moment arrays; 500 trajectories at 102001 times require about 58 GiB for those
arrays alone.

The analysis-local `post_analysis.ipynb` loads that artifact and its metadata,
averages complete trajectories with sampling covariance, constructs J moments
and modes, and retains raw moments in `mcwf_moments.generator_moments` for the
early-time trajectory bootstrap. Solve MFE
locally using the artifact's physical parameters, protocol, and saved-time
grid. Use the shared converters and plotters throughout. Plot aggregated and
group-resolved J angles and moments, bright/dark vectors, phase-2 dark traces,
and the phase-2 dark-component windowed harmonic metrics. Reuse the reference
notebook's window durations, fit-quality display criteria, and diagnostic
counts. Write MCWF and MFE steps as separate explicit calls, including their
mode construction, analysis, and plot overlays. Keep notebook analysis and
plotting cells free of loops; pass curves together and use the harmonic
quality and reference-overlay helpers for repeated work.
Keep the early-time damped fit separate and call
`compute_short_time_dark_mode_analysis(...)` with the same duration,
resampling settings, and MFE reference frequencies. Follow
`docs/instructions/post_analysis/harmonic_analysis.typ` for both analyses.
MCWF component and length plots receive matching SEMs. Follow the reference
notebook's colors and line styles: phase-2 traces distinguish MCWF and MFE by
color, while the other comparisons use dashed MFE overlays. Show only bright
and dark mode curves. Notebook Markdown cells contain titles only. Optional
figure export goes to `output/plots/<run-name>`. Include the dashed MFE dark
curve in each early-time fit panel, using the same interval and time origin.

= Invariants

- Mode construction is post-analysis and must remain independent of the solver.
- Construct each mode container from x/y/z and optional sampling uncertainty;
  `BlochVectorSeries` is authoritative for deriving lengths and SEMs.
- Do not compute modes from `nx`, `ny`, and `nz`, which
  are instantaneous unit directions. Construct inputs as
  $2 bold(J)_g / N_(J,g)$ so vector shrinkage remains visible.
- Do not silently generalize the two-group dark mode to more groups; a
  multi-group dark subspace requires an explicit basis convention.
