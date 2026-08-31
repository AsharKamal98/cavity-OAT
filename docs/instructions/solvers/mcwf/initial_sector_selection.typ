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

Support two mutually exclusive selection inputs:

```python
dN: int | None
joint_probability: float | None
linear_constraints: (A, b) | None
```

Exactly one of `dN` and `joint_probability` must be non-`None`.
`linear_constraints` is an independent optional restriction. In all cases, let

$
N = sum_(i=1)^G N_i,
quad
N_J = sum_(i=1)^G N_(J,i).
$

== Optional Linear Constraints

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

Let the constraint-eligible tuple space be

$
Omega_C = {bold(N)_J in product_(i=1)^G {0, ..., N_i}:
  A bold(N)_J = bold(b)}.
$

With no constraints, $Omega_C$ is the complete Cartesian tuple space. Reject
constraints for which $Omega_C$ is empty. These equations constrain the sector
coordinates $N_(J,i)$, not the fixed group sizes $N_i$.

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

Without `linear_constraints`, this includes every valid decomposition and
makes the retained physical total-$N_J$ subspace independent of how the same
atoms are divided into computational groups. With constraints, `dN` retains
the intersection of that complete window with $Omega_C$.

== `joint_probability`: Highest-Probability Tuples

`joint_probability` specifies a maximum retained probability $q$ relative to
the complete unconstrained binomial wavefunction, with

$
0 < q <= 1.
$

Enumerate the complete Cartesian tuple space, apply `linear_constraints`, and
compute the exact integer weight of every tuple in $Omega_C$ before selecting
any:

$
W(bold(N)_J) = product_(i=1)^G binom(N_i, N_(J,i)).
$

Define the full and constraint-eligible weights

$
Z_("full") = sum_(bold(N)_J) W(bold(N)_J) = 2^N,
quad
Z_C = sum_(bold(N)_J in Omega_C) W(bold(N)_J).
$

The eligible probability in the complete wavefunction is

$
P_C = frac(Z_C, Z_("full")).
$

If $P_C <= q$, retain every eligible tuple. Otherwise, sort the eligible
tuples by decreasing $W$ and consider complete equal-weight groups in that
order. Add a group only when

$
frac(Z_S + Z_("group"), Z_("full")) <= q,
$

where $Z_S$ is the weight already retained. Stop before the first complete
tie group that would exceed the cap. Thus the retained full-wavefunction
probability satisfies $P_S = Z_S / Z_("full") <= q$ and may be below $q$.
Never split an equal-weight group. If even the highest-weight complete tie
group exceeds $q$, reject the input because no nonempty tie-preserving support
fits under the requested cap.

Do not use a heap or best-first traversal for this method. Compute all exact
weights first, then sort them. This makes the global ordering and cutoff-tie
handling explicit and deterministic. Unlike unconstrained `dN`, this selection
acts on individual group tuples and therefore depends on the chosen group
partition.

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

Require `sector_distribution="binomial"` when `joint_probability` is used.
The square distribution is defined over complete total-$N_J$ sectors and has
no consistent meaning for a probability-ranked subset of individual tuples.

= Required Selection Summary

For either selection method, compute

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
P_C = frac(Z_C, 2^N),
quad
P_S = frac(sum_(bold(N)_J in S) W(bold(N)_J), 2^N),
quad
0 < P_S <= P_C <= 1.
$

Here $K_("ret")$ counts the group-resolved tuple blocks that will actually be
evolved, $K_C$ counts the constraint-eligible tuples, and $K_("full")$ counts
the complete Cartesian tuple space. Always print one summary during ensemble
setup. With constraints, use for example:

```text
Initial sector selection: 12/31 eligible sectors retained (38.71%; 12/961 of full space, 1.25%); full-wavefunction probability: 7.98% retained from 12.54% eligible (63.64% of eligible mass).
```

This print is required for both `dN` and `joint_probability` and must not be
gated by `verbose`. Print it once after the shared sector selection, not inside
per-trajectory state construction.

Compute the complete exact weight table once before either cutoff, then apply
the optional constraints. This provides both $Z_("full")$ and $Z_C$;
`joint_probability` additionally sorts the eligible table. For
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


def select_initial_sector_weights(
    Ni,
    *,
    dN=None,
    joint_probability=None,
    linear_constraints=None,
):
    validate exactly one of dN and joint_probability is non-None
    all_sectors = product(*(range(Ni_i + 1) for Ni_i in Ni))
    all_weights = exact_group_sector_weights(Ni, all_sectors)
    eligible_weights = constraint_eligible_weights(
        all_weights, linear_constraints
    )
    reject an empty eligible_weights dictionary

    if dN is not None:
        totals = centered_total_nj_values(sum(Ni), dN)
        retained_weights = {
            sector: weight
            for sector, weight in eligible_weights.items()
            if sum(sector) in totals
        }
        return retained_weights, eligible_weights

    validate 0 < joint_probability <= 1
    full_weight_sum = sum(all_weights.values())
    if sum(eligible_weights.values()) / full_weight_sum <= joint_probability:
        return eligible_weights, eligible_weights
    sort eligible_weights in decreasing order without a heap
    retain each complete equal-weight group only while its addition keeps
        cumulative_weight / full_weight_sum <= joint_probability
    stop before the first group that would exceed joint_probability
    reject the input if no complete tie group can be retained
    return retained_weights, eligible_weights


def compute_initial_sector_summary(Ni, selected_weights, eligible_weights):
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
    }


def print_initial_sector_summary(summary):
    print the retained count, eligible count, eligible percentage,
        retained count versus full space, full-space percentage,
        eligible and retained full-wavefunction probabilities on one line


def group_generic_initial_sector_coeffs(
    Ni,
    *,
    dN=None,
    joint_probability=None,
    linear_constraints=None,
    sector_distribution="binomial",
):
    weights, eligible_weights = select_initial_sector_weights(
        Ni,
        dN=dN,
        joint_probability=joint_probability,
        linear_constraints=linear_constraints,
    )
    summary = compute_initial_sector_summary(Ni, weights, eligible_weights)
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
        joint_probability=parameters.joint_probability,
        linear_constraints=parameters.linear_constraints,
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
`constraint_eligible_weights(...)` restricts it to $Omega_C$ before
`select_initial_sector_weights(...)` applies either cutoff.
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
- Apply `linear_constraints` before either selection method, but use the full
  unconstrained weight $2^N$ as the denominator for `joint_probability` and
  both reported probabilities.
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
- `joint_probability` never retains more than its requested fraction of the
  complete unconstrained binomial wavefunction. If the complete constrained
  support already lies below the cap, retain all of it.
- Preserve complete equal-weight groups. Reject a cap below the probability of
  the highest complete tie group rather than returning an empty state,
  exceeding the cap, or breaking the tie.
- Grouping-equivalence comparisons must use unconstrained `dN` with the same
  total $N$ and `sector_distribution`; `joint_probability` and explicit linear
  constraints do not define grouping-independent retained subspaces.
- The `dN` implementation may reject odd total $N$. Future odd-$N$ support must
  retain complete equal-probability pairs around $N/2$ rather than choose an
  asymmetric cutoff. `joint_probability` does not require even total $N$.
- Both methods enumerate $product_i (N_i + 1)$ tuples before constraints and
  selection. Validate the requested problem size before allocating the full
  weight table. Do not add multiprocessing solely for coefficient construction.
