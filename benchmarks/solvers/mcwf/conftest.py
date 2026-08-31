from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
import pytest

from common.utils.phases import default_three_phase_protocol
from parser.common import Phase, PhaseProtocol
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.mcwf import MCWFSolverParameters, SectorKey
from parser.moments import SimulationMetadata
from solvers.mcwf.sim import build_precomputed_trajectory_data
from solvers.mcwf.state_helpers import (
    build_initial_sector_state,
    group_generic_initial_sector_coeffs,
)


@dataclass(frozen=True)
class PrecomputeBenchmarkCase:
    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    integration_phases: tuple[Phase, ...]
    sector_coeffs: Mapping[SectorKey, complex]
    dt: float


@dataclass(frozen=True)
class TrajectoryBenchmarkCase(PrecomputeBenchmarkCase):
    t_eval: np.ndarray
    seed: int
    shifted_jump_operator: bool


@dataclass(frozen=True)
class PropagationBenchmarkCase:
    psi_blocks: tuple[np.ndarray, ...]
    generators: tuple[Any, ...]
    propagators: tuple[Any, ...]
    partial_step: float


@dataclass(frozen=True)
class EnsembleBenchmarkCase:
    parameters: MCWFSolverParameters
    t_eval: np.ndarray
    ntraj: int
    seed: int
    n_processes: int
    chunksize: int


@dataclass(frozen=True)
class MomentAnalysisBenchmarkCase:
    moments: TrajectoryResolvedMomentSeries
    metadata: SimulationMetadata


@pytest.fixture(scope="session")
def phase_protocol() -> PhaseProtocol:
    return default_three_phase_protocol(
        durations=(0.10, 0.15, 0.05),
        delta0=0.4,
        Omega0=4.0,
    )


@pytest.fixture(scope="session")
def integration_phases(phase_protocol: PhaseProtocol) -> tuple[Phase, ...]:
    return phase_protocol.integration_phases


@pytest.fixture(scope="session")
def jump_heavy_phase_protocol() -> PhaseProtocol:
    return default_three_phase_protocol(
        durations=(0.50, 1.00, 0.10),
        delta0=1.0,
        Omega0=16.0,
    )


@pytest.fixture(scope="session")
def precompute_case(
    integration_phases: tuple[Phase, ...],
) -> PrecomputeBenchmarkCase:
    Ni = (20, 20)
    sector_coeffs, _ = group_generic_initial_sector_coeffs(Ni, dN=2)
    return PrecomputeBenchmarkCase(
        Ni=Ni,
        omega_i=(0.9, 1.1),
        Gamma=0.1,
        integration_phases=integration_phases,
        sector_coeffs=sector_coeffs,
        dt=0.005,
    )


@pytest.fixture(scope="session")
def no_jump_phase_protocol() -> PhaseProtocol:
    return default_three_phase_protocol(
        durations=(0.10, 0.15, 0.05),
        delta0=0.4,
        Omega0=0.0,
    )


@pytest.fixture(scope="session")
def no_jump_trajectory_case(
    no_jump_phase_protocol: PhaseProtocol,
) -> TrajectoryBenchmarkCase:
    Ni = (8, 8)
    sector_coeffs, _ = group_generic_initial_sector_coeffs(Ni, dN=1)
    return TrajectoryBenchmarkCase(
        Ni=Ni,
        omega_i=(0.9, 1.1),
        Gamma=0.001,
        integration_phases=no_jump_phase_protocol.integration_phases,
        sector_coeffs=sector_coeffs,
        dt=0.005,
        t_eval=np.array([0.0, 0.30]),
        seed=1234,
        shifted_jump_operator=True,
    )


@pytest.fixture(scope="session")
def no_jump_trajectory_precomputed(no_jump_trajectory_case: TrajectoryBenchmarkCase):
    case = no_jump_trajectory_case
    return build_precomputed_trajectory_data(
        Ni=case.Ni,
        omega_i=case.omega_i,
        Gamma=case.Gamma,
        integration_phases=case.integration_phases,
        sector_coeffs=case.sector_coeffs,
        dt=case.dt,
        shifted_jump_operator=case.shifted_jump_operator,
    )


@pytest.fixture(scope="session")
def jump_heavy_trajectory_case(
    jump_heavy_phase_protocol: PhaseProtocol,
) -> TrajectoryBenchmarkCase:
    Ni = (8, 8)
    sector_coeffs, _ = group_generic_initial_sector_coeffs(Ni, dN=1)
    return TrajectoryBenchmarkCase(
        Ni=Ni,
        omega_i=(0.9, 1.1),
        Gamma=1.0,
        integration_phases=jump_heavy_phase_protocol.integration_phases,
        sector_coeffs=sector_coeffs,
        dt=0.005,
        t_eval=np.linspace(0.0, jump_heavy_phase_protocol.total_duration, 61),
        seed=1234,
        shifted_jump_operator=False,
    )


@pytest.fixture(scope="session")
def jump_heavy_trajectory_precomputed(
    jump_heavy_trajectory_case: TrajectoryBenchmarkCase,
):
    case = jump_heavy_trajectory_case
    return build_precomputed_trajectory_data(
        Ni=case.Ni,
        omega_i=case.omega_i,
        Gamma=case.Gamma,
        integration_phases=case.integration_phases,
        sector_coeffs=case.sector_coeffs,
        dt=case.dt,
        shifted_jump_operator=case.shifted_jump_operator,
    )


@pytest.fixture(scope="session")
def propagation_case(precompute_case: PrecomputeBenchmarkCase) -> PropagationBenchmarkCase:
    case = precompute_case
    precomputed = build_precomputed_trajectory_data(
        Ni=case.Ni,
        omega_i=case.omega_i,
        Gamma=case.Gamma,
        integration_phases=case.integration_phases,
        sector_coeffs=case.sector_coeffs,
        dt=case.dt,
    )
    blocks = build_initial_sector_state(case.sector_coeffs)
    phase_index = 1
    return PropagationBenchmarkCase(
        psi_blocks=tuple(blocks[key] for key in precomputed["sector_list"]),
        generators=tuple(precomputed["integration_phase_generators"][phase_index]),
        propagators=tuple(precomputed["integration_phase_propagators"][phase_index]),
        partial_step=0.5 * case.dt,
    )


@pytest.fixture(scope="session")
def ensemble_phase_protocol() -> PhaseProtocol:
    return default_three_phase_protocol(
        durations=(0.10, 0.15, 0.05),
        delta0=0.05,
        Omega0=0.5,
    )


@pytest.fixture(scope="session")
def serial_ensemble_case(
    ensemble_phase_protocol: PhaseProtocol,
) -> EnsembleBenchmarkCase:
    parameters = MCWFSolverParameters(
        Ni=(8, 8),
        omega_i=(0.9, 1.1),
        Gamma=0.2,
        phase_protocol=ensemble_phase_protocol,
        dN=1,
        dt=0.005,
    )
    return EnsembleBenchmarkCase(
        parameters=parameters,
        t_eval=np.linspace(0.0, ensemble_phase_protocol.total_duration, 31),
        ntraj=8,
        seed=2468,
        n_processes=1,
        chunksize=1,
    )


@pytest.fixture(scope="session")
def multiprocessing_ensemble_case(
    serial_ensemble_case: EnsembleBenchmarkCase,
) -> EnsembleBenchmarkCase:
    case = serial_ensemble_case
    return EnsembleBenchmarkCase(
        parameters=case.parameters,
        t_eval=case.t_eval,
        ntraj=32,
        seed=case.seed,
        n_processes=2,
        chunksize=4,
    )


@pytest.fixture(scope="session")
def moment_analysis_case(
    phase_protocol: PhaseProtocol,
) -> MomentAnalysisBenchmarkCase:
    Ni = (8, 8)
    t_eval = np.linspace(0.0, phase_protocol.total_duration, 5001)
    num_trajectories = 4
    num_groups = len(Ni)
    generator_count = 8 * num_groups
    packed_count = generator_count * (generator_count + 1) // 2
    first_moments = np.zeros(
        (num_trajectories, t_eval.size, num_groups, 8),
        dtype=float,
    )
    time_profile = np.linspace(-0.25, 0.25, t_eval.size)
    for trajectory_index in range(num_trajectories):
        first_moments[trajectory_index, :, :, 5] = (
            time_profile[:, np.newaxis] + 0.01 * trajectory_index
        )
    second_moments = np.zeros(
        (num_trajectories, t_eval.size, packed_count),
        dtype=float,
    )

    return MomentAnalysisBenchmarkCase(
        moments=TrajectoryResolvedMomentSeries(
            t=t_eval,
            first_moments=first_moments,
            second_moments=second_moments,
        ),
        metadata=SimulationMetadata(
            Ni=Ni,
            omega_i=(0.9,),
            Gamma=0.2,
            phase_protocol=phase_protocol,
        ),
    )
