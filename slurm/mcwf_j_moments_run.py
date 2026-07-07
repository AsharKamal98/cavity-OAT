from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from common.utils.parameters import default_three_phase_protocol, omega_c
from parser.mcwf import MCWFSolverParameters
from parser.moments import MomentSeries
from slurm.j_moments_io import save_j_moments_artifact
from solvers.mcwf.ensamble_sim import run_trajectory_ensemble
from solvers.mcwf.j_moments import compute_mcwf_j_moments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MCWF and save averaged J moments.")
    parser.add_argument("--n-processes", type=int, required=True)
    parser.add_argument("--ntraj", type=int, required=True)
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--array-index", type=int, required=True)
    return parser.parse_args()


def main() -> None:
    main_t0 = time.perf_counter()
    args = parse_args()

    output_dir = Path(__file__).resolve().parent / "outputs"

    # common fixed
    Gamma = 1.0
    dt = 1e-2
    num_snapshots = 100
    # param
    dN = 0
    Ni = [10, 10]
    omega_i = [0.7]
    # model and parameters
    Omega_ratio = 0.4
    Omega0 = Omega_ratio * omega_c(sum(Ni) // 2, Gamma)
    delta0 = 1.0
    # protocol phases
    phases = default_three_phase_protocol(
        T1=10.0,
        T2=10.0,
        T3=10.0,
        delta0=delta0,
        Omega0=Omega0,
    )

    total_time = float(sum(phase.duration for phase in phases))

    mcwf_moments = MomentSeries(
        num_snapshots=num_snapshots,
        total_time=total_time,
    )

    mcwf_parameters = MCWFSolverParameters(
        Ni=Ni,
        dN=dN,
        omega_i=omega_i,
        Gamma=Gamma,
        phases=phases,
        sector_distribution="binomial",
        dt=dt,
        shifted_jump_operator=True,
    )
    simulation_t0 = time.perf_counter()
    mcwf_ensemble = run_trajectory_ensemble(
        mcwf_parameters,
        t_eval=mcwf_moments.t,
        seed=1234 + args.array_index,
        ntraj=args.ntraj,
        n_processes=args.n_processes,
        verbose=True,
    )
    simulation_time = time.perf_counter() - simulation_t0
    print(f"Simulation runtime: {simulation_time:.2f} seconds.")

    j_moments_t0 = time.perf_counter()
    mcwf_moments.J = compute_mcwf_j_moments(
        mcwf_ensemble,
        n_processes=args.n_processes,
    )
    j_moments_time = time.perf_counter() - j_moments_t0
    print(f"J-moments runtime: {j_moments_time:.2f} seconds.")

    output_path = output_dir / f"{args.filename}_{args.array_index}.pkl"
    save_j_moments_artifact(mcwf_moments.J, phases, output_path)
    print(f"Saved J moments artifact to {output_path}")
    print(f"Total main runtime: {time.perf_counter() - main_t0:.2f} seconds.")


if __name__ == "__main__":
    main()
