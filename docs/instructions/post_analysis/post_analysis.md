# Post Analysis

This file does not need to follow the `agent-instruction-writer` skill format.
It should simply list functions in the `post_analysis/` package with short
summaries attached to them.

## `post_analysis/generator_moments.py`

- `compute_mean_generator_moments(moments, *, trajectory_indices=None)`
  averages trajectory-resolved first and symmetric second generator moments
  over the selected complete trajectories. The returned
  `TrajectoryResolvedMomentSeries` retains a singleton trajectory axis.
  Repeated indices are allowed for bootstrap resampling. This function does
  not construct or average covariance matrices; covariance must be formed from
  the averaged first and second moments in later analysis.

## `post_analysis/j_moments.py`

- `compute_j_moments(mean_generator_moments, *, metadata, tol=1e-12)` converts
  already-averaged first-order SU(3) generator moments into the shared
  group-resolved `JMomentSeries`. The input must retain a singleton trajectory
  axis. It uses the fixed mappings $J_x=G_6$, $J_y=-G_7$,
  $J_z=G_3/2-\sqrt{3}G_8/2$ and the corresponding linear combinations for
  $N_e$ and $N_J$, then attaches group and aggregated directions and angles
  through the `JMomentSeries` class methods. It leaves the symmetric second
  moments available for later S-moment, covariance, and squeezing analysis.
  The detailed contract is in `docs/instructions/post_analysis/j_moments.typ`.

## `post_analysis/s_moments.py`

- `compute_s_moments(mean_generator_moments, j_moments, *, tol=1e-12)`
  constructs the group-local S frame from the group-resolved J angles and
  contracts its time-dependent $3\times8$ generator-coefficient matrix with
  the already-averaged first generator moments. It returns an `SMomentSeries`
  containing raw components, the effective zero/one basis, and the lengths,
  normalized directions, and S angles attached through `SMomentSeries` class
  methods. The input generator series must retain a singleton trajectory axis
  and share its time and group axes with the J angles. Symmetric second moments
  are not transformed by this function. The detailed contract is in
  `docs/instructions/post_analysis/s_moments.typ`.

## `post_analysis/dressed_basis.py`

- `build_dressed_basis_transformation(j_moments, s_moments, *, tol=1e-12)`
  constructs the unitary, time- and group-resolved dressed basis with columns
  $(|c\rangle, |j\rangle, |s\rangle)$ in lab coordinates. It validates matching
  J/S grids and angles and returns an array shaped `(len(t), G, 3, 3)`. The
  detailed convention is in
  `docs/instructions/post_analysis/dressed_basis.typ`.

## `post_analysis/generator_rotation.py`

- `compute_gell_mann_rotation(transformation, *, tol=1e-12)` evaluates
  $R_{\alpha a}=2\operatorname{Tr}(U g_\alpha U^\dagger g_a)$ for
  $g_a=\lambda_a/2$, validates that the result is real and orthogonal, and
  returns it with shape `(len(t), G, 8, 8)`.
- `rotate_generator_moments(moments, rotation, *, tol=1e-12)` applies the
  group-local rotations to trajectory-resolved first moments and the
  block-diagonal global rotation to the full symmetric second-moment matrix.
  It preserves trajectories, cross-group blocks, and the existing packed
  upper-triangle layout. The detailed contract is in
  `docs/instructions/post_analysis/generator_rotation.typ`.

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

## `post_analysis/squeezing.py`

- `compute_squeezing(mean_dressed_generator_moments, j_moments, *,
  metadata, include_full_system=False, tol=1e-12)` is the complete
  solver-independent covariance and squeezing calculation. It consumes both
  orders of already-averaged generator moments after the caller has explicitly
  rotated them into the J/S-defined dressed frame, constructs the `4G x 4G`
  covariance of `(T1, T2, T4, T5)`, and returns a `SqueezingResult`. Its
  `group_resolved` tuple contains one
  four-observable covariance, ordered covariance eigenvalues, $N_c$, and
  generalized squeezing series per physical group. Only when
  `include_full_system=True`, it also
  sums every group-pair covariance block and constructs `full_system` for the
  unweighted sum of corresponding group fluctuation operators. The detailed
  contract is in `docs/instructions/post_analysis/squeezing.typ`.

## `post_analysis/theory_benchmarks.py`

- `phase1_ss_angles_for_nj(Nj, Omega, Gamma)` returns the phase-1
  steady-state benchmark angles `(theta_ss, phi_ss)` for one active sector.
