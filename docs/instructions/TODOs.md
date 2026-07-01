# REFACTOR TODOS

1. build_phase_jump_operator_for_sector defined in sim, but called in precompute, sim, j_moments. Documented in precompute
    1. Move function from sim to helper?
    2. Move instruction in sector operators?

2. build_t_eval_from_phases(...) is in sim.py, while MomentSeries independently rebuilds the same grid logic. Not obviously wrong, but if we want one source of truth for t_eval, this could become a small shared helper.

3. clean up _attach_mfe_residuals in `quantum_trajectoriesj_moments/.py`


# INSTRUCTION FILES

1. `sector_operators.typ`. Needs more physics, currently documented in `simulation_precompute.typ`. 

# MAYBIES

1. _interp_series(...) is duplicated in aggregator.py, squeezing.py, and theory_validation.py. This looks like a generic utility and could move to quantum_trajectories/utils.py, but it is not as central to the physics API.
