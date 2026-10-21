from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any

import numpy as np

from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.mcwf import TrajectoryEnsemble
from parser.moments import SimulationMetadata


def save_trajectory_ensemble_artifact(
    ensemble: TrajectoryEnsemble,
    simulation_metadata: SimulationMetadata,
    path: Path,
    *,
    array_index: int,
) -> None:
    """Save one array job's complete MCWF ensemble."""
    if ensemble.metadata is None:
        raise ValueError("TrajectoryEnsemble.metadata is required.")
    if ensemble.moments.num_trajectories == 0:
        raise ValueError("Cannot save an empty TrajectoryEnsemble.")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump(
            {
                "ensemble": ensemble,
                "simulation_metadata": simulation_metadata,
                "array_index": array_index,
                "ntraj": ensemble.moments.num_trajectories,
            },
            fh,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


def load_trajectory_ensemble_artifact(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        return pickle.load(fh)


def load_combined_trajectory_ensemble(
    input_dir: Path,
    filename: str,
    num_files: int,
) -> tuple[TrajectoryEnsemble, SimulationMetadata]:
    """Join raw moment chunks in array-index order without averaging."""
    if num_files < 1:
        raise ValueError("num_files must be positive.")

    paths = [input_dir / f"{filename}_{index}.pkl" for index in range(num_files)]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Missing trajectory chunk: {path}")

    simulation_metadata: SimulationMetadata | None = None
    reference_metadata = None
    trajectory_counts: list[int] = []
    seeds: list[tuple[int, ...]] = []
    for index, path in enumerate(paths):
        artifact = load_trajectory_ensemble_artifact(path)
        ensemble = artifact.get("ensemble")
        context = artifact.get("simulation_metadata")
        if not isinstance(ensemble, TrajectoryEnsemble):
            raise TypeError(f"{path} does not contain a TrajectoryEnsemble.")
        if not isinstance(context, SimulationMetadata):
            raise TypeError(f"{path} does not contain SimulationMetadata.")
        if artifact.get("array_index") != index:
            raise ValueError(f"{path} has the wrong array index.")
        if artifact.get("ntraj") != ensemble.moments.num_trajectories:
            raise ValueError(f"{path} has an inconsistent trajectory count.")
        if len(ensemble.seeds) != ensemble.moments.num_trajectories:
            raise ValueError(f"{path} has an inconsistent seed count.")
        if ensemble.metadata is None:
            raise ValueError(f"{path} has no ensemble metadata.")
        if ensemble.moments.uncertainty is not None:
            raise ValueError(f"{path} must contain raw, unaveraged trajectories.")

        metadata = ensemble.metadata
        if (
            metadata.Ni != context.Ni
            or metadata.omega_i != context.omega_groups
            or metadata.Gamma != context.Gamma
            or metadata.phase_protocol != context.phase_protocol
        ):
            raise ValueError(f"{path} has inconsistent physical metadata.")
        if not np.array_equal(metadata.t_eval, ensemble.moments.t):
            raise ValueError(f"{path} has inconsistent time grids.")

        if reference_metadata is not None:
            if context != simulation_metadata:
                raise ValueError(f"{path} has different simulation metadata.")
            for field_name in (
                "Ni", "omega_i", "Gamma", "phase_protocol",
                "shifted_jump_operator", "sectors",
                "sector_multiplicities", "sector_dimensions",
            ):
                if getattr(metadata, field_name) != getattr(
                    reference_metadata, field_name
                ):
                    raise ValueError(f"{path} has different {field_name} metadata.")
            if not np.array_equal(ensemble.moments.t, time_grid):
                raise ValueError(f"{path} has a different saved-time grid.")
        else:
            simulation_metadata = context
            reference_metadata = metadata
            time_grid = ensemble.moments.t.copy()
            first_shape = ensemble.moments.first_moments.shape[1:]
            second_shape = ensemble.moments.second_moments.shape[1:]
        if (
            ensemble.moments.first_moments.shape[1:] != first_shape
            or ensemble.moments.second_moments.shape[1:] != second_shape
        ):
            raise ValueError(f"{path} has incompatible moment axes.")
        trajectory_counts.append(ensemble.moments.num_trajectories)
        # Spawn keys restart in each independently seeded array task.
        seeds.extend((index, *seed) for seed in ensemble.seeds)
        del artifact, ensemble

    # Read each chunk twice so only one input chunk and the combined arrays
    # need to fit in memory. The first pass above validates and counts rows.
    total_trajectories = sum(trajectory_counts)
    first_moments = np.empty((total_trajectories, *first_shape), dtype=float)
    second_moments = np.empty((total_trajectories, *second_shape), dtype=float)
    start = 0
    for path, count in zip(paths, trajectory_counts):
        artifact = load_trajectory_ensemble_artifact(path)
        moments = artifact["ensemble"].moments
        stop = start + count
        first_moments[start:stop] = moments.first_moments
        second_moments[start:stop] = moments.second_moments
        start = stop
        del artifact, moments

    combined = TrajectoryEnsemble(
        moments=TrajectoryResolvedMomentSeries(
            t=time_grid,
            first_moments=first_moments,
            second_moments=second_moments,
        ),
        seeds=seeds,
        metadata=reference_metadata,
    )
    assert simulation_metadata is not None
    return combined, simulation_metadata
