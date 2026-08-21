from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np

from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.mcwf import TrajectoryEnsemble, TrajectoryEnsembleMetadata
from parser.moments import MomentSeries, SimulationMetadata
from post_analysis import (
    build_dressed_basis_transformation,
    compute_gell_mann_rotation,
    compute_j_moments,
    compute_mean_generator_moments,
    compute_s_moments,
    compute_squeezing,
    rotate_generator_moments,
)
from slurm.j_moments_io import save_mcwf_moments_artifact
from slurm.trajectory_ensemble_io import load_trajectory_ensemble_artifact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine MCWF ensembles and compute squeezing moments."
    )
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--num-files", type=int, default=4)
    return parser.parse_args()


def _validate_matching_metadata(
    metadata: TrajectoryEnsembleMetadata,
    reference: TrajectoryEnsembleMetadata,
) -> None:
    fields = (
        "Ni",
        "omega_i",
        "Gamma",
        "phase_protocol",
        "shifted_jump_operator",
        "sectors",
        "sector_multiplicities",
        "sector_dimensions",
    )
    if any(
        getattr(metadata, field) != getattr(reference, field)
        for field in fields
    ):
        raise ValueError("All ensemble files must share the same metadata.")
    if not np.allclose(metadata.t_eval, reference.t_eval, atol=1e-12, rtol=0.0):
        raise ValueError("All ensemble files must share the same t_eval grid.")


def _load_combined_ensemble(
    input_dir: Path,
    filename: str,
    num_files: int,
) -> tuple[TrajectoryEnsemble, SimulationMetadata]:
    if num_files <= 0:
        raise ValueError("num_files must be positive.")

    first_moment_chunks = []
    second_moment_chunks = []
    seeds = []
    simulation_metadata_ref = None
    ensemble_metadata_ref = None

    for array_index in range(num_files):
        path = input_dir / f"{filename}_{array_index}.pkl"
        artifact = load_trajectory_ensemble_artifact(path)
        required_fields = (
            "ensemble",
            "simulation_metadata",
            "array_index",
            "ntraj",
        )
        missing_fields = [field for field in required_fields if field not in artifact]
        if missing_fields:
            raise ValueError(
                f"{path} is missing: " + ", ".join(missing_fields)
            )

        ensemble = artifact["ensemble"]
        simulation_metadata = artifact["simulation_metadata"]
        if not isinstance(ensemble, TrajectoryEnsemble):
            raise TypeError(f"{path} does not contain a TrajectoryEnsemble.")
        if not isinstance(simulation_metadata, SimulationMetadata):
            raise TypeError(f"{path} does not contain SimulationMetadata.")
        if artifact["array_index"] != array_index:
            raise ValueError(f"{path} has the wrong array index.")
        if artifact["ntraj"] != ensemble.moments.num_trajectories:
            raise ValueError(f"{path} has inconsistent trajectory counts.")
        if len(ensemble.seeds) != ensemble.moments.num_trajectories:
            raise ValueError(f"{path} has inconsistent seed counts.")
        if ensemble.metadata is None:
            raise ValueError(f"{path} has no ensemble metadata.")

        if simulation_metadata_ref is None:
            simulation_metadata_ref = simulation_metadata
            ensemble_metadata_ref = ensemble.metadata
        else:
            if simulation_metadata != simulation_metadata_ref:
                raise ValueError(
                    "All ensemble files must share the same simulation metadata."
                )
            _validate_matching_metadata(ensemble.metadata, ensemble_metadata_ref)

        first_moment_chunks.append(ensemble.moments.first_moments)
        second_moment_chunks.append(ensemble.moments.second_moments)
        seeds.extend((array_index, *seed) for seed in ensemble.seeds)

    if simulation_metadata_ref is None or ensemble_metadata_ref is None:
        raise ValueError("No trajectory ensemble files were loaded.")

    return (
        TrajectoryEnsemble(
            moments=TrajectoryResolvedMomentSeries(
                t=np.asarray(ensemble_metadata_ref.t_eval, dtype=float).copy(),
                first_moments=np.concatenate(first_moment_chunks, axis=0),
                second_moments=np.concatenate(second_moment_chunks, axis=0),
            ),
            seeds=seeds,
            metadata=ensemble_metadata_ref,
        ),
        simulation_metadata_ref,
    )


def _compute_moments(
    ensemble: TrajectoryEnsemble,
    simulation_metadata: SimulationMetadata,
) -> MomentSeries:
    moments = MomentSeries(
        t=ensemble.moments.t,
        metadata=simulation_metadata,
    )
    moments.generator_moments = ensemble.moments
    moments.mean_generator_moments = compute_mean_generator_moments(
        moments.generator_moments,
    )
    moments.J = compute_j_moments(
        moments.mean_generator_moments,
        metadata=moments.metadata,
    )
    moments.S = compute_s_moments(
        moments.mean_generator_moments,
        moments.J,
    )
    transformation = build_dressed_basis_transformation(moments.J, moments.S)
    rotation = compute_gell_mann_rotation(transformation)
    mean_dressed_generator_moments = rotate_generator_moments(
        moments.mean_generator_moments,
        rotation,
    )
    moments.squeezing = compute_squeezing(
        mean_dressed_generator_moments,
        moments.J,
        metadata=moments.metadata,
        include_full_system=False,
    )
    return moments


def main() -> None:
    args = parse_args()
    ensemble, simulation_metadata = _load_combined_ensemble(
        args.input_dir,
        args.filename,
        args.num_files,
    )
    print(f"Loaded {ensemble.moments.num_trajectories} trajectories.")

    analysis_t0 = time.perf_counter()
    moments = _compute_moments(ensemble, simulation_metadata)
    print(f"Post-analysis runtime: {time.perf_counter() - analysis_t0:.2f} seconds.")

    output_path = args.output_dir / f"{args.filename}.pkl"
    save_mcwf_moments_artifact(
        moments,
        output_path,
        ntraj=ensemble.moments.num_trajectories,
    )
    print(f"Saved MCWF moments to {output_path}")


if __name__ == "__main__":
    main()
