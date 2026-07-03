# REFACTOR TODOS

## General
1. build_phase_jump_operator_for_sector defined in sim, but called in precompute, sim, j_moments. Documented in precompute
    1. Move function from sim to helper?
    2. Move instruction in sector operators?

2. build_t_eval_from_phases(...) is in sim.py, while MomentSeries independently rebuilds the same grid logic. Not obviously wrong, but if we want one source of truth for t_eval, this could become a small shared helper.

3. clean up _attach_mfe_residuals in `quantum_trajectoriesj_moments/.py`

4. Go through legacy functions in helpers in `quantum_trajectories`

5. Rename 
    `quantum_trajectories` -> `mcwf`
    `quantum_trajectories/qutip` -> `qutip_fixed_nj`

6. Move mfe_plotting to common

7. refactor sim functions in `quantum_trajectories`

8. MomentSeries Parameters make no sense. 

## Homoheneous - Inhomogeneous split to group-split
Want to stop treating simulations and post-soimulation analysis differently depending on if we have homogeneous or inhomogeneous setup. Split should be based on number of groups, independent of what effective couplings omega_i they have.

1. Input to simulation should N_i = [N1,N2,...NG], omega_i = [omega1,omega2,...,omegaG]. Remove N input.
2. precompute constructs operators over Hilbert space split into G smaller Hilbert spaces. This can be done later, for now hard-coded two-group split works fine.
3. simulator itself should be easily extendable, since it already works with generic tuple sectors for inhomogeneous case. Can be fixed later.
4. Compute J moment series should be changed so group fields are default, and filled even for G=1. The non-group fields become extra fields contining averaged moments.
5. Plotting accepts lists of series (groups), or single series (averages). No need to distinguish between homogeneous or inhomogeneous runs.

## MFE Residuals
1. Move mfe_residuals to post-analysis. 
2. Create class method which uses residual util and attatches it to class.

# INSTRUCTION FILES

1. `sector_operators.typ`. Needs more physics, currently documented in `simulation_precompute.typ`. 

# MAYBIES

1. _interp_series(...) is duplicated in aggregator.py, squeezing.py, and theory_validation.py. This looks like a generic utility and could move to quantum_trajectories/utils.py, but it is not as central to the physics API.
