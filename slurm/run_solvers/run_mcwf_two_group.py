from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from common.utils.parameters import scaled_N_Gamma
from common.utils.phases import default_three_phase_protocol
from parser.moments import SimulationMetadata
from slurm.run_solvers.run_mcwf_common import (
    add_common_arguments,
    parse_float_list,
    parse_int_list,
    run_mcwf_case,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run two-group MCWF and save averaged J moments.")
    return add_common_arguments(parser).parse_args()


def main() -> None:
    main_t0 = time.perf_counter()
    args = parse_args()

    output_dir = args.output_dir

    Gamma = args.Gamma
    dt = args.dt
    num_snapshots = args.num_snapshots
    dN = args.dN
    Ni = parse_int_list(args.Ni)
    omega_i = parse_float_list(args.omega_i)
    N = sum(Ni)
    Omega0 = scaled_N_Gamma(args.Omega_factor, N, Gamma)
    delta0 = scaled_N_Gamma(args.delta_factor, N, Gamma)
    phase_protocol = default_three_phase_protocol(
        durations=(args.T1, args.T2, args.T3),
        delta0=delta0,
        Omega0=Omega0,
    )
    metadata = SimulationMetadata(
        Ni=tuple(Ni),
        omega_i=tuple(omega_i),
        Gamma=Gamma,
        phase_protocol=phase_protocol,
    )

    run_mcwf_case(
        label="two_group",
        output_path=output_dir / f"{args.filename}_{args.array_index}.pkl",
        metadata=metadata,
        dN=dN,
        dt=dt,
        num_snapshots=num_snapshots,
        seed=1234 + args.array_index,
        ntraj=args.ntraj,
        n_processes=args.n_processes,
    )

    print(f"Total two_group main runtime: {time.perf_counter() - main_t0:.2f} seconds.")


if __name__ == "__main__":
    main()
