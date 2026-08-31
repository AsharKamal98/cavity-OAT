#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Group-Generic Initial Sector Selection and State]]

= Purpose

This file specifies the target group-generic selection of initial MCWF
strong-symmetry sectors, their coefficients, and their internal state blocks.
Use it when editing the group-generic helpers in
`solvers/mcwf/state_helpers.py` or the initialization flow in
`solvers/mcwf/ensamble_sim.py` and `solvers/mcwf/sim.py`.
The binomial initial-state convention follows the coherent product-state
expansion in `docs/theory/main.tex`.

Represent the full initial wavefunction as a direct sum of sector blocks. The
MCWF simulation should evolve one internal array per sector in the
group-resolved active-manifold excitation basis.

= Selection Contract

== Public Inputs

Support independent sector filters:

```python
dN: int | None
max_omega_over_omega_c: float | None
linear_constraints: (A, b) | None
joint_probability: float | None
```

All four filters are independent and optional. Execute a filter exactly when
its input is non-`None`; if all are `None`, retain the complete tuple support.
Apply the construction in distinct ordered stages:

1. enumerate the complete tuple space and compute exact weights;
2. apply the complete-total-$N_J$ `dN` filter, when supplied;
3. apply the `max_omega_over_omega_c` drive-validity filter, when supplied;
4. apply the $N_J$-ratio equations in `linear_constraints`, when supplied;
5. apply `joint_probability`, when supplied, as the final support filter;
6. normalize the retained coefficients and build the sector blocks.

In all cases, let

$
N = sum_(i=1)^G N_i,
quad
N_J = sum_(i=1)^G N_(J,i).
$

== `dN`: Complete Total-$N_J$ Window

`dN` specifies the half-width of the retained total-$N_J$ window around
$N/2$.
For even $N$, require

$
0 <= d N <= N/2.
$

The retained totals are

$
T = {N/2-d N, ..., N/2+d N},
quad
M = |T| = 2 d N + 1.
$

Thus `dN=2` retains five complete total-$N_J$ sectors. The final number of
group-resolved tuple blocks is generally larger than $M$.

Represent a group-resolved sector by the tuple

$
bold(N)_J = (N_(J,1), ..., N_(J,G)).
$

For every retained $N_J in T$, include every constraint-eligible tuple
satisfying

$
sum_(i=1)^G N_(J,i) = N_J,
quad
0 <= N_(J,i) <= N_i.
$

At this stage, include every valid decomposition. This makes the `dN` filter
independent of how the same atoms are divided into computational groups. The
following ratio filter may then restrict this support.

If `dN` is `None`, define $Omega_D$ as the complete Cartesian tuple space. If
it is supplied, let $Omega_D$ be the tuples in the complete space whose total
$N_J$ lies in $T$.

== Maximum Drive-Ratio Filter

`max_omega_over_omega_c` specifies the strict upper bound $r_(max)$ on the
drive ratio, with

$
0 < r_(max) <= 1,
quad
frac(Omega_(max), Omega_c(N_J)) < r_(max),
quad
Omega_c(N_J) = frac(Gamma N_J, 2),
$

where $Omega_(max)$ is the largest absolute drive over every integration
phase. Compute

$
N_(J,min)
=
floor(frac(2 Omega_(max), Gamma r_(max))) + 1,
$

then retain every tuple in $Omega_D$ satisfying

$
sum_i N_(J,i) >= N_(J,min).
$

This removes complete total-$N_J$ shells and therefore preserves every group
decomposition that survived `dN`. If `max_omega_over_omega_c` is `None`, skip
this filter. Call the existing drive-ratio check after selection as a final
validation using the same $Omega_(max)$ and configured limit.

Let $Omega_R$ denote the support after this stage, with $Omega_R=Omega_D$ when
the filter is omitted.

== Optional $N_J$-Ratio Constraints

Use `linear_constraints=(A, b)` to retain only group-resolved tuples satisfying

$
A bold(N)_J = bold(b).
$

`A` must be an integer matrix with $G$ columns, and `b` must be an integer
vector with one entry per row of `A`. Apply the equations exactly, without a
floating-point tolerance. For example,

```python
# NJ1 = 2 NJ2
linear_constraints=(
    ((1, -2),),
    (0,),
)

# NJ1 = 2 NJ2 + 3 NJ3
linear_constraints=(
    ((1, -2, -3),),
    (0,),
)
```

The second example defines a plane in three-dimensional sector space. A true
line in $G$ dimensions generally requires $G-1$ independent equations.

Apply these equations to $Omega_R$, after the optional `dN` and maximum
drive-ratio filters. Let

$
Omega_C = {bold(N)_J in Omega_R: A bold(N)_J = bold(b)}.
$

With no constraints, $Omega_C=Omega_R$. Reject constraints for which $Omega_C$
is empty. These equations constrain the sector coordinates $N_(J,i)$, not the
fixed group sizes $N_i$.

== `joint_probability`: Final Tuple-Wise Probability Filter

`joint_probability` specifies a maximum retained probability $q$ relative to
the complete unconstrained binomial wavefunction, with

$
0 < q <= 1.
$

Compute the exact integer weight of every tuple before applying any filters:

$
W(bold(N)_J) = product_(i=1)^G binom(N_i, N_(J,i)).
$

Define the full weight

$
Z_("full") = sum_(bold(N)_J) W(bold(N)_J) = 2^N.
$

Let $Omega_F$ be the candidate tuple set remaining after `dN`, the maximum
drive-ratio filter, and $N_J$-ratio constraints. The probability before the
final filter is

$
P_F = frac(sum_(bold(N)_J in Omega_F) W(bold(N)_J), Z_("full")).
$

If $P_F <= q$, retain every candidate tuple. Otherwise, sort the candidate
tuples by decreasing $W$ and consider complete equal-weight groups in that
order. Add a group only when

$
frac(Z_S + Z_("group"), Z_("full")) <= q,
$

where $Z_S$ is the weight already retained. Stop before the first complete
tie group that would exceed the cap. Thus the retained full-wavefunction
probability satisfies $P_S = Z_S / Z_("full") <= q$ and may be below $q$.
Never split an exact equal-weight tie group. If even the highest-weight
complete tie group exceeds $q$, reject the input because no nonempty
tie-preserving support fits under the requested cap.

This filter always acts on individual group-resolved tuples
$bold(N)_J=(N_(J,1), ..., N_(J,G))$. It does *not* retain or remove complete
sets of tuples sharing one total $N_J$. Consequently, it may retain one tuple
with a given total $N_J$ while removing another tuple with that same total.
Equal-weight ties are handled together only to make an exactly degenerate
cutoff deterministic; an equal-weight tie group is not a total-$N_J$ sector.

Do not use a heap or best-first traversal for this filter. Compute all exact
weights first, then sort the remaining tuple candidates. This makes the global
ordering and cutoff-tie handling explicit and deterministic. Unlike
unconstrained `dN`, this filter depends on the chosen group partition.

= Exact Sector Coefficients

== Binomial Initial State

For the initial product state with equal single-atom amplitudes in $|u chevron.r$
and $|d chevron.r$, the exact group-tuple probability is

$
P(bold(N)_J)
=
frac(W(bold(N)_J), 2^N).
$

If $S$ is the set of all retained tuples and

$
Z_S = sum_(bold(N)_J in S) W(bold(N)_J),
$

the normalized sector coefficient is

$
c_(bold(N)_J)
= sqrt(frac(W(bold(N)_J), Z_S)).
$

Precompute each group's exact integer binomial table and reuse it across
tuples. Do not replace these weights with a Gaussian approximation: the exact
one-dimensional tables are inexpensive, remain valid in the tails, and avoid
confusing probabilities with wavefunction amplitudes.

For unconstrained `dN`, the retained full-wavefunction probability can also be
computed without
reference to the group split,

$
P_S = sum_(N_J in T) frac(binom(N, N_J), 2^N),
$

using the product-binomial convolution identity. Record this value as a
truncation diagnostic when the result container provides such metadata.

== Square Distribution Compatibility

For `dN`, if `sector_distribution="square"` remains supported, preserve its
current meaning: every retained total-$N_J$ sector has equal total probability,
while its group decompositions follow their conditional binomial
probabilities. For $M = |T|$,

$
|c_(bold(N)_J)|^2
=
frac(1, M)
frac(
  product_(i=1)^G binom(N_i, N_(J,i)),
  binom(N, N_J)
).
$

Require `sector_distribution="binomial"` when `joint_probability` is used,
including when it follows `dN`.
The square distribution is defined over complete total-$N_J$ sectors and has
no consistent meaning for a probability-ranked subset of individual tuples.

= Required Selection Summary

For the complete staged selection, compute

$
K_("full") = product_(i=1)^G (N_i + 1),
quad
K_C = |Omega_C|,
quad
K_("ret") = |S|,
quad
p_("sectors") = 100 frac(K_("ret"), K_C),
$

the eligible and retained probabilities in the complete unconstrained
wavefunction before coefficient renormalization,

$
Z_C = sum_(bold(N)_J in Omega_C) W(bold(N)_J),
quad
P_C = frac(Z_C, 2^N),
quad
P_S = frac(sum_(bold(N)_J in S) W(bold(N)_J), 2^N),
quad
0 < P_S <= P_C <= 1.
$

Here $K_("ret")$ counts the group-resolved tuple blocks that will actually be
evolved, $K_C$ counts the tuples remaining after `dN`, maximum drive-ratio,
and $N_J$-ratio filtering, and $K_("full")$ counts the complete Cartesian tuple
space. Keep the probability fields in `InitialSectorSelectionSummary` as
diagnostics even though the standard console table reports sector counts.

For every applied filter, store an `InitialSectorFilterSummary` containing its
name, input sector count, and remaining sector count. Always print the initial
total followed by one aligned row per applied filter. For example:

```text
Initial sector construction: 961 sectors.

Filter                       Sectors remaining   Retained from input   Retained from initial
--------------------------------------------------------------------------------------------
dN filter                                  541                56.30%                  56.30%
Omega/Omega_c filter                       400                73.94%                  41.62%
NJ-ratio filter                              8                 2.00%                   0.83%
Joint-probability filter                     3                37.50%                   0.31%
```

`Retained from input` uses the count passed into that filter. `Retained from
initial` uses $K_("full")$. Omitted filters produce no row. If every filter is
omitted, print only the initial total. This output must not be gated by
`verbose`; print it once after shared sector selection, not inside
per-trajectory state construction.

Compute the complete exact weight table once before filtering, then apply the
filters in their defined order. This provides both $Z_("full")$ and $Z_C$;
`joint_probability` additionally sorts only the tuples remaining after all
preceding filters. For
`sector_distribution="square"`, $P_S$ is still the full-wavefunction binomial
probability contained in the selected support, not the renormalized
square-distribution mass.

= Sector Blocks

Use tuple keys internally for every group count, including `(NJ,)` for one
group. A compatibility adapter may accept historical scalar one-group keys at
the public boundary.

For tuple $bold(N)_J$, the flattened internal excitation block has size

$
D_(bold(N)_J) = product_(i=1)^G (N_(J,i) + 1).
$

The default initial internal state is the all-active-atoms-down basis state:
allocate a zero block of this size and place $c_(bold(N)_J)$ at flattened index
zero.

The resulting direct-sum state must satisfy

$
sum_(bold(N)_J in S) norm(psi_(bold(N)_J))^2 = 1.
$

= Method in Pseudo-code

```python
def centered_total_nj_values(total_N, dN):
    validate even total_N
    validate 0 <= dN <= total_N // 2
    return tuple(range(total_N // 2 - dN,
                       total_N // 2 + dN + 1))


def exact_group_sector_weights(Ni, sectors):
    binomial_tables = precompute exact binomial(Ni, NJi) integers
    return {
        sector: product(table[NJi] for table, NJi in zip(
            binomial_tables, sector
        ))
        for sector in sectors
    }


def constraint_eligible_weights(all_weights, linear_constraints):
    if linear_constraints is None:
        return all_weights
    A, b = validate integer constraint dimensions
    return {
        sector: weight
        for sector, weight in all_weights.items()
        if A @ sector == b
    }


def filter_sector_weights_by_total_nj(candidate_weights, total_N, dN):
    totals = centered_total_nj_values(total_N, dN)
    return {
        sector: weight
        for sector, weight in candidate_weights.items()
        if sum(sector) in totals
    }


def minimum_total_nj_for_omega_ratio(
    max_abs_omega,
    Gamma,
    max_omega_over_omega_c,
):
    validate max_abs_omega >= 0, Gamma > 0
    validate 0 < max_omega_over_omega_c <= 1
    threshold = 2 * max_abs_omega / (
        Gamma * max_omega_over_omega_c
    )
    return floor(threshold) + 1


def filter_sector_weights_by_omega_ratio(
    candidate_weights,
    max_abs_omega,
    Gamma,
    max_omega_over_omega_c,
):
    minimum_total_nj = minimum_total_nj_for_omega_ratio(
        max_abs_omega,
        Gamma,
        max_omega_over_omega_c,
    )
    return {
        sector: weight
        for sector, weight in candidate_weights.items()
        if sum(sector) >= minimum_total_nj
    }


def filter_sector_weights_by_joint_probability(
    candidate_weights,
    joint_probability,
    full_weight_sum,
):
    validate 0 < joint_probability <= 1
    if sum(candidate_weights.values()) / full_weight_sum <= joint_probability:
        return candidate_weights
    sort individual candidate tuples in decreasing weight order without a heap
    retain each complete exact equal-weight tie group only while its addition
        keeps cumulative_weight / full_weight_sum <= joint_probability
    stop before the first tie group that would exceed joint_probability
    reject the input if no complete tie group can be retained
    return retained_weights


def _select_initial_sector_weights_with_steps(
    Ni,
    *,
    dN=None,
    max_omega_over_omega_c=None,
    max_abs_omega=None,
    Gamma=None,
    linear_constraints=None,
    joint_probability=None,
):
    all_sectors = product(*(range(Ni_i + 1) for Ni_i in Ni))
    all_weights = exact_group_sector_weights(Ni, all_sectors)
    retained_weights = all_weights
    filter_steps = []

    if dN is not None:
        input_count = len(retained_weights)
        retained_weights = filter_sector_weights_by_total_nj(
            retained_weights,
            total_N=sum(Ni),
            dN=dN,
        )
        append ("dN filter", input_count, len(retained_weights)) to filter_steps

    if max_omega_over_omega_c is not None:
        input_count = len(retained_weights)
        retained_weights = filter_sector_weights_by_omega_ratio(
            retained_weights,
            max_abs_omega,
            Gamma,
            max_omega_over_omega_c,
        )
        append ("Omega/Omega_c filter", input_count,
                len(retained_weights)) to filter_steps

    if linear_constraints is not None:
        input_count = len(retained_weights)
        retained_weights = constraint_eligible_weights(
            retained_weights, linear_constraints
        )
        append ("NJ-ratio filter", input_count,
                len(retained_weights)) to filter_steps

    pre_probability_weights = retained_weights

    if joint_probability is not None:
        input_count = len(retained_weights)
        retained_weights = filter_sector_weights_by_joint_probability(
            retained_weights,
            joint_probability,
            full_weight_sum=sum(all_weights.values()),
        )
        append ("Joint-probability filter", input_count,
                len(retained_weights)) to filter_steps

    reject an empty retained_weights dictionary
    return retained_weights, pre_probability_weights, filter_steps


def compute_initial_sector_summary(
    Ni,
    selected_weights,
    eligible_weights,
    filter_steps,
):
    total_sector_count = product(Ni_i + 1 for Ni_i in Ni)
    eligible_sector_count = len(eligible_weights)
    retained_sector_count = len(selected_weights)
    return {
        "retained_sector_count": retained_sector_count,
        "eligible_sector_count": eligible_sector_count,
        "total_sector_count": total_sector_count,
        "retained_percentage_of_eligible": (
            100 * retained_sector_count / eligible_sector_count
        ),
        "retained_percentage_of_full_space": (
            100 * retained_sector_count / total_sector_count
        ),
        "eligible_joint_probability": (
            sum(eligible_weights.values()) / 2 ** sum(Ni)
        ),
        "retained_joint_probability": (
            sum(selected_weights.values()) / 2 ** sum(Ni)
        ),
        "filter_steps": filter_steps,
    }


def print_initial_sector_summary(summary):
    print the initial total sector count
    print an aligned header
    for step in summary.filter_steps:
        print step name, remaining count,
            100 * remaining / input_count,
            100 * remaining / initial_count


def group_generic_initial_sector_coeffs(
    Ni,
    *,
    dN=None,
    max_omega_over_omega_c=None,
    max_abs_omega=None,
    Gamma=None,
    linear_constraints=None,
    joint_probability=None,
    sector_distribution="binomial",
):
    weights, pre_probability_weights, filter_steps = (
        _select_initial_sector_weights_with_steps(
        Ni,
        dN=dN,
        max_omega_over_omega_c=max_omega_over_omega_c,
        max_abs_omega=max_abs_omega,
        Gamma=Gamma,
        linear_constraints=linear_constraints,
        joint_probability=joint_probability,
        )
    )
    summary = compute_initial_sector_summary(
        Ni, weights, pre_probability_weights, filter_steps
    )
    coefficients = coefficients_from_selected_weights(
        weights,
        sector_distribution=sector_distribution,
        selection_uses_dN=dN is not None,
    )
    return normalize_sector_coefficients(coefficients), summary


def build_initial_sector_state(sector_coeffs):
    blocks = {}
    for sector, coefficient in sector_coeffs.items():
        block = zeros(product(NJi + 1 for NJi in sector), complex128)
        block[0] = coefficient
        blocks[sector] = block
    return blocks


def _simulate_single_trajectory(..., sector_coeffs, precomputed):
    initial_sector_blocks = build_initial_sector_state(
        sector_coeffs
    )
    evolve initial_sector_blocks using precomputed data
    return TrajectoryResult(...)


def run_trajectory_ensemble(parameters, *, t_eval, ntraj, ...):
    sector_coeffs, sector_summary = group_generic_initial_sector_coeffs(
        parameters.Ni,
        dN=parameters.dN,
        max_omega_over_omega_c=parameters.max_omega_over_omega_c,
        max_abs_omega=max(abs(phase.omega) for phase in integration_phases),
        Gamma=parameters.Gamma,
        linear_constraints=parameters.linear_constraints,
        joint_probability=parameters.joint_probability,
        sector_distribution=parameters.sector_distribution,
    )
    print_initial_sector_summary(sector_summary)  # always print once
    check_initial_sector_omega_ratio(...)
    precomputed = build_precomputed_trajectory_data(
        ..., sector_coeffs=sector_coeffs
    )
    simulate every trajectory with _simulate_single_trajectory(
        ..., sector_coeffs=sector_coeffs, precomputed=precomputed
    )
    return TrajectoryEnsemble(...)
```

`centered_total_nj_values(...)` defines the `dN` window.
`exact_group_sector_weights(...)` computes the complete exact tuple table, and
`filter_sector_weights_by_total_nj(...)` applies the first optional
complete-sector window. `minimum_total_nj_for_omega_ratio(...)` derives the
strict lower bound used by `filter_sector_weights_by_omega_ratio(...)`, which
acts immediately afterward. `constraint_eligible_weights(...)` then applies
the optional exact $N_J$-ratio equations, while
`filter_sector_weights_by_joint_probability(...)` applies the final cap
to individual remaining tuples.
`_select_initial_sector_weights_with_steps(...)` owns this ordered filter
pipeline and records one `InitialSectorFilterSummary` per applied filter. The
public `select_initial_sector_weights(...)` may preserve its existing two-item
return by discarding only these reporting rows.
`group_generic_initial_sector_coeffs(...)` converts the retained weights into
normalized amplitudes and returns their summary. The summary is computed by
`compute_initial_sector_summary(...)` and printed once by
`print_initial_sector_summary(...)`. `build_initial_sector_state(...)` then
converts the amplitudes into the direct-sum blocks consumed by
`_simulate_single_trajectory(...)`.
Initial-sector setup belongs inside `run_trajectory_ensemble(...)`, not in a
notebook: the ensemble entry point selects and validates the shared sector
coefficients once, constructs reusable precomputed data, and gives every
trajectory the same sector representation.

= Implementation Scope

- Use the group-generic coefficient flow above instead of dispatching between
  separate one- and two-group initial-state helpers.
- Do not reinterpret `dN` as a number of individual group tuples; it always
  selects a complete total-$N_J$ window.
- Apply `dN` first, `max_omega_over_omega_c` second, `linear_constraints`
  third, and `joint_probability` last. Execute each filter only when its input
  is non-`None`. Use the full unconstrained weight
  $2^N$ as the denominator for `joint_probability` and both reported
  probabilities.
- Keep this implementation under `solvers/mcwf/`; no other current solver
  consumes the direct-sum sector-wavefunction representation.
- A group-generic initializer does not by itself enable $G > 2$ propagation.
  Sector operators, multiplicities, and precomputation must also support
  arbitrary-length tuple keys according to their own MCWF instructions.

= Invariants and Edge Cases

- Initial construction is deterministic for fixed inputs.
- The printed selection summary is unconditional and appears exactly once per
  ensemble setup, independent of trajectory count and `verbose`.
- Sector coefficients satisfy $sum_(bold(N)_J in S) |c_(bold(N)_J)|^2 = 1$.
- The represented sector keys and sector metadata are identical for every
  trajectory in an ensemble.
- The internal excitation basis and flattened block shape are fixed by the
  sector key and must not change during a trajectory.
- Without linear constraints, every selected total-$N_J$ value includes all
  valid group decompositions.
- `max_omega_over_omega_c` removes complete total-$N_J$ shells below the
  strict lower bound computed from the largest absolute integration-phase
  drive. The final drive-ratio check must use the same limit.
- `joint_probability` never retains more than its requested fraction of the
  complete unconstrained binomial wavefunction. If the support remaining after
  all preceding filters already lies below the cap, retain all of it.
- `joint_probability` filters individual group-resolved tuples. It never
  changes into a complete-total-$N_J$ filter, even when `dN` was applied first.
  Tuples sharing one total $N_J$ may therefore receive different outcomes.
- Preserve complete equal-weight groups. Reject a cap below the probability of
  the highest complete tie group rather than returning an empty state,
  exceeding the cap, or breaking the tie.
- Grouping-equivalence comparisons must use unconstrained `dN` with the same
  total $N$ and `sector_distribution`; the tuple-wise `joint_probability`
  filter and explicit linear constraints do not define grouping-independent
  retained subspaces.
- The `dN` implementation may reject odd total $N$. Future odd-$N$ support must
  retain complete equal-probability pairs around $N/2$ rather than choose an
  asymmetric cutoff. `joint_probability` does not require even total $N$.
- The staged constructor enumerates $product_i (N_i + 1)$ tuples before
  constraints and filtration. Validate the requested problem size before
  allocating the full weight table. Do not add multiprocessing solely for
  coefficient construction.
