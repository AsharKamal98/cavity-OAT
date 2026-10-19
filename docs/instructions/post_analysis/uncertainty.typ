#set page(paper: "a4", margin: 18mm)
#set text(size: 10pt)
#set par(first-line-indent: 0pt, spacing: 0.65em)
#set heading(numbering: "1.")
#show raw: set text(size: 9pt)

#align(center)[#text(size: 16pt, weight: "bold")[Moment Uncertainty]]

= Purpose

Define the shared sampling-covariance and standard-deviation flow in
`common/utils/uncertainty.py`. Existing imports through
`post_analysis/uncertainty.py` and `post_analysis` remain available as exports.
Observable modules and attachment class methods supply their own
transformations and attach the results; their instruction files define the
physics and output layout.

= General Flow

At a fixed time, the covariance matrix contains each variable's variance on
the diagonal and the covariances between variables elsewhere. For two variables,

$
C_x = mat(
  s_1^2, C_(12);
  C_(12), s_2^2;
),
quad
C_(i j) = frac(1, R-1) sum_(r=1)^R
  (x_(r,i)-overline(x)_i)(x_(r,j)-overline(x)_j).
$

Here $R$ is the number of independent trajectories, $overline(x)_i$ is the
trajectory average, and $s_i$ is the trajectory standard deviation.

For a fixed affine relation $bold(y)=A bold(x)+bold(b)$, $A$ is the
transformation matrix: row $i$ contains the coefficients defining output
$y_i$. Propagate covariance and extract standard deviations and SEMs using

$
C_y = A C_x A^T,
quad s_y = sqrt("diag"(C_y)),
quad "SEM"_y = s_y / sqrt(R).
$

The square root acts entrywise; the constant offset $bold(b)$ does not affect
covariance. For a nonlinear relation $bold(y)=f(bold(x))$, use the Jacobian
(the matrix of first derivatives) evaluated at the mean input:

$
A_(i j) = frac(partial f_i, partial x_j)
quad "evaluated at" quad bold(x)=overline(bold(x)).
$

The same covariance formula is then a first-order approximation. Evaluate
the Jacobian separately at each saved time using that time's averaged inputs.

```python
stats = compute_sample_covariance(samples, mean=sample_mean)
if stats is not None:
    covariance = propagate_covariance(stats.sample_covariance, A)
    std = std_from_covariance(covariance)
    sem = std / sqrt(stats.n_samples)
```

= Two-Variable Example

Suppose $x_1$ and $x_2$ have standard deviations $s_1=2$ and $s_2=1$, with
covariance $C_(12)=1$. Their covariance matrix is

$
C_x = mat(4, 1; 1, 1).
$

Define $y_1=x_1+x_2$ and $y_2=x_1-x_2$. Their means obey the same relations:
$overline(y)_1=overline(x)_1+overline(x)_2$ and
$overline(y)_2=overline(x)_1-overline(x)_2$. The transformation is

$
vec(y_1, y_2) =
underbrace(mat(1, 1; 1, -1), A) vec(x_1, x_2).
$

These coefficients are also the Jacobian entries; they are constant, so the
propagation is exact. Multiplying the matrices gives

$
C_y = A C_x A^T
= mat(1, 1; 1, -1) mat(4, 1; 1, 1) mat(1, 1; 1, -1)
= mat(7, 3; 3, 3).
$

The diagonal entries are the new variances, and the off-diagonal entry is
the new covariance. Thus

$
s_(y_1) = sqrt(7) approx 2.65,
quad s_(y_2) = sqrt(3) approx 1.73,
quad "Cov"(y_1,y_2)=3.
$

Equivalently, expanding the squared deviations gives
$s_(y_1)^2=4+1+2 times 1=7$ and
$s_(y_2)^2=4+1-2 times 1=3$: positive covariance increases the uncertainty of
the sum and reduces that of the difference. The two original standard
deviations alone would not supply this covariance contribution.

For the averages over $R$ independent trajectories,

$
"SEM"_(y_1)=sqrt(frac(7,R)),
quad "SEM"_(y_2)=sqrt(frac(3,R)).
$

For example, $R=100$ gives SEMs of approximately $0.265$ and $0.173$.

= Helper Contracts

- `compute_sample_covariance(...)` computes the covariance matrix across
  trajectories and returns it with the trajectory count.
- `propagate_covariance(...)` computes the new variables' covariance matrix
  as $A C A^T$ using the supplied transformation matrix or Jacobian.
- `std_from_covariance(...)` computes each variable's standard deviation by
  taking the square root of its covariance matrix's diagonal entry.

= Storage and Responsibilities

`SamplingUncertainty` in `parser/uncertainty.py` contains only
`sample_covariance` and `n_samples`. Generator averaging automatically
attaches it as `mean_generator_moments.uncertainty`, using one full
`(time, 8G, 8G)` matrix with cross-group correlations. Averaging and selection
conventions are documented in `docs/instructions/parser.typ` and
`docs/instructions/post_analysis/post_analysis.md`.
With only one trajectory, sample covariance cannot be estimated and
`compute_sample_covariance(...)` returns `None`.

Observable extraction functions and attachment class methods should compute
their central values and delegate associated standard deviations to these
helpers. Missing input
uncertainty leaves output uncertainty unset. Retain full covariance for
subsequent conversions and aggregation; separate standard deviations do not
encode correlations. The caller converts propagated sample standard
deviations to SEMs by dividing by $sqrt(R)$ and attaches the results.
Evaluate derivatives only at valid points and mark undefined derived SEMs
with `NaN`; the stored source covariance remains finite.
These statistics describe trajectory sampling, not quantum fluctuations.
