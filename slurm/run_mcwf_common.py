from __future__ import annotations

import argparse
import time
from pathlib import Path

from parser.mcwf import MCWFSolverParameters
from parser.moments import MomentSeries
from slurm.j_moments_io import save_j_moments_artifact
from solvers.mcwf.ensamble_sim import run_trajectory_ensemble
from solvers.mcwf.j_moments import compute_mcwf_j_moments


def add_common_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_omega_i: bool = True,
) -> argparse.ArgumentParser:
    parser.add_argument("--n-processes", type=int, required=True)
    parser.add_argument("--ntraj", type=int, required=True)
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--array-index", type=int, required=True)
    parser.add_argument("--Gamma", type=float, required=True)
    parser.add_argument("--Ni", type=str, required=True)
    parser.add_argument("--dN", type=int, required=True)
    if include_omega_i:
        parser.add_argument("--omega-i", type=str, required=True)
    parser.add_argument("--Omega-factor", dest="Omega_factor", type=float, required=True)
    parser.add_argument("--delta-factor", dest="delta_factor", type=float, required=True)
    parser.add_argument("--T1", type=float, required=True)
    parser.add_argument("--T2", type=float, required=True)
    parser.add_argument("--T3", type=float, required=True)
    return parser


def parse_float_list(value: str) -> list[float]:
    if not value.strip():
        return []
    return [float(item.strip()) for item in value.split(",")]


def parse_int_list(value: str) -> list[int]:
    if not value.strip():
        return []
    return [int(item.strip()) for item in value.split(",")]


def run_mcwf_case(
    *,
    label: str,
    output_path: Path,
    Ni: list[int],
    omega_i: list[float],
    dN: int,
    Gamma: float,
    phases,
    dt: float,
    num_snapshots: int,
    seed: int,
    ntraj: int,
    n_processes: int,
) -> None:
    print(f"Starting {label} MCWF run with Ni={Ni}, omega_i={omega_i}.")
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
        seed=seed,
        ntraj=ntraj,
        n_processes=n_processes,
        verbose=True,
    )
    simulation_time = time.perf_counter() - simulation_t0
    print(f"{label} simulation runtime: {simulation_time:.2f} seconds.")

    j_moments_t0 = time.perf_counter()
    mcwf_moments.J = compute_mcwf_j_moments(
        mcwf_ensemble,
        n_processes=n_processes,
    )
    j_moments_time = time.perf_counter() - j_moments_t0
    print(f"{label} J-moments runtime: {j_moments_time:.2f} seconds.")

    save_j_moments_artifact(mcwf_moments.J, phases, output_path)
    print(f"Saved {label} J moments artifact to {output_path}")
