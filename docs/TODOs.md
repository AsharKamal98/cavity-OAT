# REFACTOR TODOS


## Refactors

1. MCWF Homoheneous - Inhomogeneous split to group-split
    Want to stop treating simulations and post-soimulation analysis differently depending on if we have homogeneous or inhomogeneous setup. Split should be based on number of groups G, independent of what effective couplings omega_i they have.
    1. DONE. Input to simulation should N_i = [N1,N2,...NG], omega_i = [omega1,omega2,...,omega_G-1]. Remove N input.
    2. precompute constructs operators over Hilbert space split into G smaller Hilbert spaces. This can be done later, for now hard-coded two-group split works fine.
    3. simulator itself should be easily extendable, since it already works with generic tuple sectors for inhomogeneous case. Can be fixed later.
    4. Compute J moment series should be changed so group fields are default, and filled even for G=1. The non-group fields become extra fields contining averaged moments.
    5. Plotting accepts lists of series (groups), or single series (averages). No need to distinguish between homogeneous or inhomogeneous runs.


2. Qutip mcsolve to group-split
    1. DONE. Ni and omega_i lists as input.
    2. Change build_qutip_two_group_fixed_nj_model_from_phases output to lists.
    3. Should print as now, but change where and how parameters accessed.
    4. Make qutip j moment computation use Ni and omega_i lists. 

3. Qutip mesolve to group-split
    1. Add Ni and omega_i dependencies.
    2. Add print statement like mcsolve.

4. j-moment computations. 
    Dependent on above.
    1. x_groups, nx_groups, N_groups etc. -> x, nx, N should become primary fields and are filled even for G=1, while x, nx etc. -> x_avg, nx_avg become secondary average files.
    2. j-moment computations in mcwf,mfe and qutip should output above.

5. `common/plotting/j_spin.py` 
    1. option 1: mid-level function. Extracts groups-resolved moments and plots in a lopp. option 2: low-level function. moment series given and plots.
    2. has option for plotting averages, only required if function mid-level 


## New stuff

1. Add back `post-analysis/squeezing.py` and `common/plotting/squeezing.py`. Exist in legacy.

2. time-dependent Omega during phase 1.

3. Create a purity post-analysis as dephasing diagnostic. 

4. Add single-group mean-field residuals to `post-analysis/mfe.py` for diagnostics and plotting.

5. Add single-group mean-field solutions to `solvers/mfe.py`.

6. Can `solvers/qutio_singel_sector.py`support arbritrary number of groups G for simulating G groups using Qutip?


## Clean-up

1. refactor sim functions in `solvers/mcwf`

2. go through legacy functions in helpers in `solvers/mcwf`

3. build_phase_jump_operator_for_sector defined in sim, but called in precompute, sim, j_moments. Documented in precompute
    1. Move function from sim to helper?
    2. Move instruction in sector operators?

4. build_t_eval_from_phases(...) is in sim.py, while MomentSeries independently rebuilds the same grid logic. Not obviously wrong, but if we want one source of truth for t_eval, this could become a small shared helper.

5. Unify QuTiP `tlist` with the repository-wide `t_eval` naming so all solver backends use the same public time-grid convention.

6. rename qutip num_points to num_snapshots?

7. Remove MomentSeries Parameters?


## INSTRUCTION FILES

1. clean up
    1. EVERYTHING

2. add back from legacy
    1. squeezing.typ
    2. dephasing_diagnostics.typ?

3. `sector_operators.typ`. Needs more physics, currently documented in `simulation_precompute.typ`. 

