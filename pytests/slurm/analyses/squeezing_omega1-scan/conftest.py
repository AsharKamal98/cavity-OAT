from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from common.utils.phases import default_three_phase_protocol
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.mcwf import MCWFSolverParameters, TrajectoryEnsemble
from parser.moments import SimulationMetadata
from slurm.trajectory_ensemble_io import save_trajectory_ensemble_artifact
from solvers.mcwf.ensamble_sim import run_trajectory_ensemble


@dataclass(frozen=True)
class SqueezingPostAnalysisCase:
    ensemble: TrajectoryEnsemble
    simulation_metadata: SimulationMetadata
    input_dir: Path
    filename: str
    chunk_sizes: tuple[int, ...]


@pytest.fixture(scope="session")
def squeezing_post_analysis_case(
    tmp_path_factory: pytest.TempPathFactory,
) -> SqueezingPostAnalysisCase:
    phase_protocol = default_three_phase_protocol(
        durations=(0.1, 0.1, 0.1),
        delta0=0.04,
        Omega0=0.2,
    )
    simulation_metadata = SimulationMetadata(
        Ni=(1, 1),
        omega_i=(1.0,),
        Gamma=1.0,
        phase_protocol=phase_protocol,
    )
    t_eval = np.linspace(0.0, phase_protocol.total_duration, 4)
    parameters = MCWFSolverParameters(
        Ni=simulation_metadata.Ni,
        omega_i=simulation_metadata.omega_groups,
        Gamma=simulation_metadata.Gamma,
        phase_protocol=phase_protocol,
        dN=0,
        dt=0.05,
    )
    simulated = run_trajectory_ensemble(
        parameters,
        t_eval=t_eval,
        seed=1234,
        ntraj=3,
        n_processes=1,
    )
    ensemble = simulated

    input_dir = tmp_path_factory.mktemp("squeezing-post-analysis")
    filename = "tiny_ensemble"
    chunk_sizes = (1, 2)
    chunk_start = 0
    for array_index, chunk_size in enumerate(chunk_sizes):
        chunk_end = chunk_start + chunk_size
        chunk = TrajectoryEnsemble(
            moments=TrajectoryResolvedMomentSeries(
                t=ensemble.moments.t.copy(),
                first_moments=ensemble.moments.first_moments[
                    chunk_start:chunk_end
                ].copy(),
                second_moments=ensemble.moments.second_moments[
                    chunk_start:chunk_end
                ].copy(),
            ),
            seeds=ensemble.seeds[chunk_start:chunk_end],
            metadata=ensemble.metadata,
        )
        save_trajectory_ensemble_artifact(
            chunk,
            simulation_metadata,
            input_dir / f"{filename}_{array_index}.pkl",
            array_index=array_index,
        )
        chunk_start = chunk_end

    return SqueezingPostAnalysisCase(
        ensemble=ensemble,
        simulation_metadata=simulation_metadata,
        input_dir=input_dir,
        filename=filename,
        chunk_sizes=chunk_sizes,
    )
