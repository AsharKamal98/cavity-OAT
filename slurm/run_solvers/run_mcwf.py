from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from common.utils.parameters import mcwf_dt_from_scales, scaled_N_Gamma
from common.utils.phases import default_three_phase_protocol
from parser.mcwf import MCWFSolverParameters
from parser.moments import MomentSeries, SimulationMetadata
from slurm.j_moments_io import save_j_moments_artifact
from slurm.trajectory_ensemble_io import save_trajectory_ensemble_artifact
from solvers.mcwf.ensamble_sim import run_trajectory_ensemble
from solvers.mcwf.j_moments import compute_mcwf_j_moments


def _parse_float_list(value: str) -> tuple[float, ...]:
    if not value.strip():
        return ()
    return tuple(float(item.strip()) for item in value.split(","))


def _parse_int_list(value: str) -> tuple[int, ...]:
    if not value.strip():
        return ()
    return tuple(int(item.strip()) for item in value.split(","))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one- or two-group MCWF and save its output."
    )
    parser.add_argument("--n-processes", type=int, required=True)
    parser.add_argument("--ntraj", type=int, required=True)
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--array-index", type=int, required=True)
    parser.add_argument("--Gamma", type=float, required=True)
    parser.add_argument("--Ni", type=_parse_int_list, required=True)
    parser.add_argument("--dN", type=int, required=True)
    parser.add_argument("--dt", type=float, default=None)
    parser.add_argument("--num-snapshots", type=int, required=True)
    parser.add_argument(
        "--omega-i",
        type=_parse_float_list,
        default=(),
        help="Comma-separated independent group couplings; omit for one group.",
    )
    parser.add_argument(
        "--Omega-factor",
        dest="Omega_factor",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--delta-factor",
        dest="delta_factor",
        type=float,
        required=True,
    )
    parser.add_argument("--T1", type=float, required=True)
    parser.add_argument("--T2", type=float, required=True)
    parser.add_argument("--T3", type=float, required=True)
    parser.add_argument(
        "--save-ensemble",
        action="store_true",
        help="Save the full ensemble instead of computing J moments.",
    )
    return parser.parse_args()


def _build_metadata(args: argparse.Namespace) -> SimulationMetadata:
    if len(args.Ni) not in (1, 2):
        raise ValueError("MCWF runs currently require one or two groups in Ni.")

    N = sum(args.Ni)
    phase_protocol = default_three_phase_protocol(
        durations=(args.T1, args.T2, args.T3),
        delta0=scaled_N_Gamma(args.delta_factor, N, args.Gamma),
        Omega0=scaled_N_Gamma(args.Omega_factor, N, args.Gamma),
    )
    return SimulationMetadata(
        Ni=args.Ni,
        omega_i=args.omega_i,
        Gamma=args.Gamma,
        phase_protocol=phase_protocol,
    )


def _resolve_dt(
    requested_dt: float | None,
    metadata: SimulationMetadata,
) -> float:
    if requested_dt is not None:
        return requested_dt

    integration_phases = metadata.phase_protocol.integration_phases
    max_omega = max(abs(phase.omega) for phase in integration_phases)
    max_delta = max(abs(phase.delta) for phase in integration_phases)
    return mcwf_dt_from_scales(
        max_omega,
        max_delta,
        sum(metadata.Ni),
        metadata.Gamma,
    )


def _run_mcwf(
    args: argparse.Namespace,
    metadata: SimulationMetadata,
) -> None:
    label = "single_group" if len(metadata.Ni) == 1 else "two_group"
    output_path = args.output_dir / f"{args.filename}_{args.array_index}.pkl"
    dt = _resolve_dt(args.dt, metadata)

    print(
        f"Starting {label} MCWF run with "
        f"Ni={metadata.Ni}, omega_i={metadata.omega_i}."
    )
    moments = MomentSeries(
        num_snapshots=args.num_snapshots,
        metadata=metadata,
    )
    parameters = MCWFSolverParameters(
        Ni=metadata.Ni,
        omega_i=metadata.omega_groups,
        Gamma=metadata.Gamma,
        phase_protocol=metadata.phase_protocol,
        dN=args.dN,
        sector_distribution="binomial",
        dt=dt,
        shifted_jump_operator=True,
    )

    simulation_t0 = time.perf_counter()
    ensemble = run_trajectory_ensemble(
        parameters,
        t_eval=moments.t,
        seed=1234 + args.array_index,
        ntraj=args.ntraj,
        n_processes=args.n_processes,
        verbose=True,
    )
    simulation_time = time.perf_counter() - simulation_t0
    print(f"{label} simulation runtime: {simulation_time:.2f} seconds.")

    if args.save_ensemble:
        save_trajectory_ensemble_artifact(
            ensemble,
            metadata,
            output_path,
            array_index=args.array_index,
        )
        print(f"Saved {label} trajectory ensemble to {output_path}")
        return

    j_moments_t0 = time.perf_counter()
    moments.J = compute_mcwf_j_moments(
        ensemble,
        n_processes=args.n_processes,
    )
    j_moments_time = time.perf_counter() - j_moments_t0
    print(f"{label} J-moments runtime: {j_moments_time:.2f} seconds.")

    save_j_moments_artifact(
        moments.J,
        metadata.phase_protocol,
        output_path,
    )
    print(f"Saved {label} moments artifact to {output_path}")


def main() -> None:
    main_t0 = time.perf_counter()
    args = parse_args()
    metadata = _build_metadata(args)
    _run_mcwf(args, metadata)
    print(f"Total MCWF runtime: {time.perf_counter() - main_t0:.2f} seconds.")


if __name__ == "__main__":
    main()
