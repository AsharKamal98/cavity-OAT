#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Two-Group J-Vector Modes]
]

= Purpose

This file defines the one- and two-group common, contrast, drive-bright, and drive-dark
vectors computed in `post_analysis/j_modes.py`. Use these derived modes for
post-simulation comparison and plotting without modifying raw `moments.J` data.

= Definitions

The inputs are the active-population-normalized group vectors

$
bold(s)_g(t) = frac(2 bold(J)_g(t), N_(J,g)(t)).
$

For two groups with active populations $N_(J,1)$ and $N_(J,2)$, define

$
bold(s)_"common" = frac(N_(J,1) bold(s)_1 + N_(J,2) bold(s)_2, N_(J,1) + N_(J,2)),
quad
bold(s)_"contrast" = frac(bold(s)_1 - bold(s)_2, 2).
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

For one group, $bold(s)_"common"=bold(s)_"bright"=bold(s)_1$ and the
contrast and dark vectors are zero.

= Method in Pseudo-code

```python
def compute_j_modes(t, x_components, y_components, z_components,
                    *, populations, omega_groups) -> JModeSeries:
    require matching x/y/z component series for one or two groups
    common = population-weighted average of the group vectors
    contrast = half difference of the group vectors
    bright, dark = coupling-weighted orthogonal combinations
    return JModeSeries(t, common, contrast, bright, dark)
```

`compute_j_modes(...)` consumes separate normalized x/y/z group series,
constructs the two vectors internally, and returns all four modes as one
`JModeSeries`.

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
J-moment representation.

= Data Requirements

- One shared time grid and matching normalized x/y/z component series for one
  or two groups. Each component input may be one series or a collection of
  group series.
- Two active-population series and two completed couplings `omega_groups`.

= Invariants

- Mode construction is post-analysis and must remain independent of the solver.
- Construct each mode container from x/y/z only; `BlochVectorSeries` is
  authoritative for deriving `length` and `xy_length`.
- Do not compute modes from `nx_groups`, `ny_groups`, and `nz_groups`, which
  are instantaneous unit directions. Construct inputs as
  $2 bold(J)_g / N_(J,g)$ so vector shrinkage remains visible.
- Do not silently generalize the two-group dark mode to more groups; a
  multi-group dark subspace requires an explicit basis convention.
