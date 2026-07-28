# REFACTOR TODOS


## Refactors

1. MCWF generalized to G groups
    Want to stop treating simulations and post-soimulation analysis differently depending on if we have homogeneous or inhomogeneous setup. Split should be based on number of groups G, independent of what effective couplings omega_i they have.
    1. DONE. Input to simulation should N_i = [N1,N2,...NG], omega_i = [omega1,omega2,...,omega_G-1]. Remove N input.
    2. precompute constructs operators over Hilbert space split into G smaller Hilbert spaces. This can be done later, for now hard-coded two-group split works fine.
    3. simulator itself should be easily extendable, since it already works with generic tuple sectors for inhomogeneous case. Can be fixed later.

2. j-moment computations. 
    Dependent on above.
    1. DONE. Group fields such as `x`, `nx`, and `N_j` are the primary
       tuple-valued fields, including for $G=1$; full-system fields live in
       the nested `aggregated` container.
    2. DONE. J-moment computations in MCWF, MFE, and QuTiP output the
       group-first field layout above.


## New stuff

1. Add back `post-analysis/squeezing.py` and `common/plotting/squeezing.py`. Exist in legacy.

3. Create a purity post-analysis as dephasing diagnostic. 

4. Add single-group mean-field residuals to `post-analysis/mfe_residuals.py` for diagnostics and plotting.


## Clean-up

1. refactor sim functions in `solvers/mcwf`

2. go through legacy functions in helpers in `solvers/mcwf`

3. build_phase_jump_operator_for_sector defined in sim, but called in precompute, sim, j_moments. Documented in precompute
    1. Move function from sim to helper?
    2. Move instruction in sector operators?


5. rename qutip num_points to num_snapshots?

6. Remove MomentSeries Parameters?

7. DONE. Investigated precomputed dyadic bisection propagators; see
   `docs/investigations/dyadic_bisection_propagators.md`.


## INSTRUCTION FILES


1. add back from legacy
    1. squeezing.typ

2. `sector_operators.typ`. Needs more physics, currently documented in `simulation_precompute.typ`.
