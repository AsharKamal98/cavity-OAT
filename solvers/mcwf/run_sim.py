from common.utils_parameters import (
    check_initial_sector_omega_ratio,
    validated_mcwf_dt,
)

import sys
import time
from solvers.mcwf.state_helpers import (
    centered_sector_initial_coeffs,
    centered_group_resolved_sector_initial_coeffs,
)
from solvers.mcwf.ensamble_sim import run_trajectory_ensemble

def run_h_sim(
    N=20,
    dN=0,
    dt=1e-2,
    Gamma=1.0,
    phases=None,
    shifted_jump_operator=True,
    ntraj=100,
    num_snapshots=100,
    seed=1234,
    n_processes=10,
):

    # sector coeff.
    homogeneous_sector_coeffs = centered_sector_initial_coeffs(
            N,
            dN=dN,
            sector_distribution="binomial",
        )
    homogeneous_ratio_check = check_initial_sector_omega_ratio(
            homogeneous_sector_coeffs,
            Omega=max(abs(phase.omega) for phase in phases),
            Gamma=Gamma,
        )
    if not homogeneous_ratio_check["is_valid"]:
        sys.exit(
            "Omega/Omega_c check not valid for homogeneous run: "
            f"Omega={homogeneous_ratio_check['omega']}, Omega_c={homogeneous_ratio_check['omega_c']}, "
            f"smallest Nj={homogeneous_ratio_check['min_nj']}, ratio={homogeneous_ratio_check['ratio']}"
        )
        return

    # simulation
    t0 = time.perf_counter()
    homogeneous_ensemble = run_trajectory_ensemble(
            N=N,
            Gamma=Gamma,
            phases=phases,
            sector_coeffs=homogeneous_sector_coeffs,
            dt=dt,
            num_snapshots=num_snapshots,
            seed=seed,
            ntraj=ntraj,
            shifted_jump_operator=shifted_jump_operator,
            n_processes=n_processes,
            chunksize=1,
            verbose=True,
        )
    simulation_time = time.perf_counter() - t0
    print("homogeneous simulation runtime:", simulation_time)

    return homogeneous_ensemble


def run_inh_sim(
    N=20,
    dN=0,
    N1=10,
    omega_1=0.7,
    dt=1e-2,
    Gamma=1.0,
    phases=None,
    shifted_jump_operator=True,
    ntraj=100,
    num_snapshots=100,
    seed=1234,
    n_processes=10,
):

    N2 = N-N1

    # sector coeff.
    inhomogeneous_sector_coeffs = centered_group_resolved_sector_initial_coeffs(
        N,
        dN=dN,
        N1=N1,
        N2=N2,
        sector_distribution="binomial",
    )
    inhomogeneous_ratio_check = check_initial_sector_omega_ratio(
            inhomogeneous_sector_coeffs,
            Omega=max(abs(phase.omega) for phase in phases),
            Gamma=Gamma,
        )
    if not inhomogeneous_ratio_check["is_valid"]:
        sys.exit(
            "Omega/Omega_c check not valid for inhomogeneous run: "
            f"Omega={inhomogeneous_ratio_check['omega']}, Omega_c={inhomogeneous_ratio_check['omega_c']}, "
            f"smallest Nj={inhomogeneous_ratio_check['min_nj']}, ratio={inhomogeneous_ratio_check['ratio']}"
        )

    t0 = time.perf_counter()
    inhomogeneous_ensemble = run_trajectory_ensemble(
        N=N,
        Gamma=Gamma,
        phases=phases,
        sector_coeffs=inhomogeneous_sector_coeffs,
        dt=dt,
        num_snapshots=num_snapshots,
        seed=seed,
        ntraj=ntraj,
        shifted_jump_operator=shifted_jump_operator,
        omega_1=omega_1,
        N1=N1,
        N2=N2,
        n_processes=n_processes,
        chunksize=1,
        verbose=True,
    )
    simulation_time = time.perf_counter() - t0
    print("inhomogeneous simulation runtime:", simulation_time)

    return inhomogeneous_ensemble
