from __future__ import annotations

import argparse
from pathlib import Path

from common.plotting import plot_bloch_angles, plot_spin_components
from common.utils.parameters import default_three_phase_protocol, omega_c
from parser.moments import MomentSeries
from slurm.j_moments_io import load_j_moments_pickle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot combined MCWF J moments.")
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--input-dir", type=Path, default=Path("slurm/outputs"))
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    Gamma = 1.0
    Ni = [10, 10]
    Omega_ratio = 0.4
    Omega0 = Omega_ratio * omega_c(sum(Ni) // 2, Gamma)
    delta0 = 1.0
    phases = default_three_phase_protocol(
        T1=10.0,
        T2=10.0,
        T3=10.0,
        delta0=delta0,
        Omega0=Omega0,
    )

    j_moments = load_j_moments_pickle(args.input_dir / f"{args.filename}.pkl")
    moments = MomentSeries(t=j_moments.t, J=j_moments)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_bloch_angles(
        moments.J,
        phases=phases,
        label="MCWF",
        colour_index=0,
        linestyle="-",
        output_path=args.output_dir / f"{args.filename}_angles.png",
    )
    plot_spin_components(
        moments.J,
        normalized=False,
        phases=phases,
        label="MCWF",
        colour_index=0,
        linestyle="-",
        output_path=args.output_dir / f"{args.filename}_spin_components.png",
    )


if __name__ == "__main__":
    main()
