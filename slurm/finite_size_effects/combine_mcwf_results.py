from __future__ import annotations

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np

from parser.j_moments import JMomentSeries
from slurm.j_moments_io import load_j_moments_artifact, save_j_moments_artifact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine array-job J-moment JSON files.")
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--num-files", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    output_dir = Path(__file__).resolve().parent / "outputs"
    artifacts = [
        load_j_moments_artifact(output_dir / f"{args.filename}_{array_index}.pkl")
        for array_index in range(args.num_files)
    ]
    samples = [artifact["J"] for artifact in artifacts]
    phases_ref = artifacts[0]["phases"]

    t_ref = np.asarray(samples[0].t, dtype=float)
    phase_ref = np.asarray(samples[0].phase_index, dtype=int)
    for artifact, sample in zip(artifacts, samples):
        if artifact["phases"] != phases_ref:
            raise ValueError("All J-moment files must share the same phases.")
        if not np.allclose(sample.t, t_ref, atol=1e-12, rtol=0.0):
            raise ValueError("All J-moment files must share the same t grid.")
        if not np.array_equal(sample.phase_index, phase_ref):
            raise ValueError("All J-moment files must share the same phase_index grid.")

    def mean_series(field_name: str):
        return np.mean(
            np.asarray([getattr(sample, field_name) for sample in samples], dtype=float),
            axis=0,
        )

    def mean_group_series(field_name: str):
        groups_ref = getattr(samples[0], field_name)
        if groups_ref is None:
            return None
        group_count = len(groups_ref)
        return tuple(
            np.mean(
                np.asarray([getattr(sample, field_name)[g] for sample in samples], dtype=float),
                axis=0,
            )
            for g in range(group_count)
        )

    combined = JMomentSeries(
        t=t_ref,
        phase_index=phase_ref,
        x=mean_series("x"),
        y=mean_series("y"),
        z=mean_series("z"),
        x_groups=mean_group_series("x_groups"),
        y_groups=mean_group_series("y_groups"),
        z_groups=mean_group_series("z_groups"),
    )
    JMomentSeries.attatch_norm_spin_components_from_spin_components(combined)
    JMomentSeries.attatch_angles_from_norm_spin_components(combined)

    output_path = output_dir / f"{args.filename}.pkl"
    save_j_moments_artifact(combined, phases_ref, output_path)
    print(f"Saved combined J moments artifact to {output_path}")


if __name__ == "__main__":
    main()
