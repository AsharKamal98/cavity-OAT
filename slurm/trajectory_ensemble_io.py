from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any

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
    if not ensemble.trajectories:
        raise ValueError("Cannot save an empty TrajectoryEnsemble.")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump(
            {
                "ensemble": ensemble,
                "simulation_metadata": simulation_metadata,
                "array_index": array_index,
                "ntraj": len(ensemble.trajectories),
            },
            fh,
        )


def load_trajectory_ensemble_artifact(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        return pickle.load(fh)
