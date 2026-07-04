# REFACTOR TODOS

## General

1. MomentSeries Parameters make no sense. 

## Clean-up

1. refactor sim functions in `solvers/mcwf`

2. go through legacy functions in helpers in `solvers/mcwf`

3. build_phase_jump_operator_for_sector defined in sim, but called in precompute, sim, j_moments. Documented in precompute
    1. Move function from sim to helper?
    2. Move instruction in sector operators?

4. build_t_eval_from_phases(...) is in sim.py, while MomentSeries independently rebuilds the same grid logic. Not obviously wrong, but if we want one source of truth for t_eval, this could become a small shared helper.

## Homoheneous - Inhomogeneous split to group-split
Want to stop treating simulations and post-soimulation analysis differently depending on if we have homogeneous or inhomogeneous setup. Split should be based on number of groups, independent of what effective couplings omega_i they have.

1. Input to simulation should N_i = [N1,N2,...NG], omega_i = [omega1,omega2,...,omegaG]. Remove N input.
2. precompute constructs operators over Hilbert space split into G smaller Hilbert spaces. This can be done later, for now hard-coded two-group split works fine.
3. simulator itself should be easily extendable, since it already works with generic tuple sectors for inhomogeneous case. Can be fixed later.
4. Compute J moment series should be changed so group fields are default, and filled even for G=1. The non-group fields become extra fields contining averaged moments.
5. Plotting accepts lists of series (groups), or single series (averages). No need to distinguish between homogeneous or inhomogeneous runs.


## INSTRUCTION FILES

1. `sector_operators.typ`. Needs more physics, currently documented in `simulation_precompute.typ`. 
