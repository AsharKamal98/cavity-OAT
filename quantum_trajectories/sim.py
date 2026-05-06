from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm, expm_multiply

from quantum_trajectories.operator_helpers import (
    build_sector_ops,
    sector_multiplicity,
)
from quantum_trajectories.parser import (
    Array,
    Phase,
    SectorOperators,
    TrajectoryResult,
    TrajectorySnapshot,
)
from quantum_trajectories.state_helpers import (
    build_initial_sector_state,
    total_norm2,
)


# -----------------------------------------------------------------------------
# Effective generators and propagation
# -----------------------------------------------------------------------------


def heff_for_sector(ops: SectorOperators, omega: float, delta: float, gamma: float) -> csc_matrix:
    """
    Regular H_delta and regular jump operator from the paper:
        H_delta = Omega J_x - delta N_e,
        l = J_-,
        H_eff = H_delta - i (gamma/2) J_+ J_-.

    In MCQT, between jumps the state evolves with the non-Hermitian effective Hamiltonian
        H_eff = H_sys - i/2 sum_l l^dagger l
    For the homogeneous collective decay, the only jump operator is J_-, so the sum reduces to J_+ J_-.
        H_delta - i/2 gamma J_+ J_-,
    """
    H = omega * ops.J_x - delta * ops.N_e
    Heff = H - 0.5j * gamma * ops.JpJm
    return Heff.tocsc()


def jump_for_sector(ops: SectorOperators, psi: Array) -> Array:
    return ops.J_minus.dot(psi)


def blocks_list_to_dict(sector_list: Sequence[int], psi_blocks: Sequence[Array]) -> Dict[int, Array]:
    return {Nj: psi.copy() for Nj, psi in zip(sector_list, psi_blocks)}


def total_norm2_list(psi_blocks: Sequence[Array]) -> float:
    return float(sum(np.vdot(psi, psi).real for psi in psi_blocks))


def renormalize_blocks(blocks: Dict[int, Array]) -> Dict[int, Array]:
    nrm = np.sqrt(total_norm2(blocks))
    if nrm == 0.0:
        raise RuntimeError("Wavefunction has zero norm.")
    return {Nj: psi / nrm for Nj, psi in blocks.items()}


def renormalize_psi_blocks(psi_blocks: Sequence[Array]) -> List[Array]:
    nrm = np.sqrt(total_norm2_list(psi_blocks))
    if nrm == 0.0:
        raise RuntimeError("Wavefunction has zero norm.")
    return [psi / nrm for psi in psi_blocks]


def propagate_blocks(
    psi_blocks: Sequence[Array],
    generators_list: Sequence[csc_matrix],
    dt: float,
) -> List[Array]:
    """
    Propagate each sector block forward in time by dt under its effective generator using expm_multiply.
    """
    out: List[Array] = []
    for psi, generator in zip(psi_blocks, generators_list):
        if psi.size == 0:
            out.append(psi.copy())
            continue
        out.append(expm_multiply((-1j * generator) * dt, psi))
    return out


def propagate_blocks_with_propagators(
    psi_blocks: Sequence[Array],
    propagators_list: Sequence[csc_matrix],
) -> List[Array]:
    """
    Propagate each sector block with a precomputed full-step propagator.
    """
    out: List[Array] = []
    for psi, propagator in zip(psi_blocks, propagators_list):
        if psi.size == 0:
            out.append(psi.copy())
            continue
        out.append(propagator @ psi)
    return out


def apply_jump(psi_blocks: Sequence[Array], ops_list: Sequence[SectorOperators]) -> List[Array]:
    jumped = [jump_for_sector(ops, psi) for ops, psi in zip(ops_list, psi_blocks)]
    if total_norm2_list(jumped) == 0.0:
        # No further jumps possible; keep the pre-jump state instead of crashing.
        return [psi.copy() for psi in psi_blocks]
    return renormalize_psi_blocks(jumped)


def build_precomputed_trajectory_data(
    N: int,
    gamma: float,
    phases: Sequence[Phase],
    sector_coeffs: Mapping[int, complex],
    dt: float,
) -> Dict[str, Any]:
    sector_list = sorted(sector_coeffs.keys())
    ops_list = [build_sector_ops(Nj) for Nj in sector_list]
    multiplicities = {Nj: sector_multiplicity(N, Nj) for Nj in sector_list}
    dims = {Nj: Nj + 1 for Nj in sector_list}
    phase_generators = [
        [heff_for_sector(ops, phase.omega, phase.delta, gamma) for ops in ops_list]
        for phase in phases
    ]
    phase_propagators = [
        [expm((-1j * generator) * dt).tocsc() for generator in generators_list]
        for generators_list in phase_generators
    ]

    return {
        "sector_list": sector_list,
        "ops_list": ops_list,
        "multiplicities": multiplicities,
        "dims": dims,
        "phase_generators": phase_generators,
        "phase_propagators": phase_propagators,
    }


# -----------------------------------------------------------------------------
# Main Monte Carlo trajectory engine
# -----------------------------------------------------------------------------


def simulate_single_trajectory(
    N: int,
    gamma: float,
    phases: Sequence[Phase],
    sector_coeffs: Mapping[int, complex],
    *,
    internal_sector_states: Optional[Mapping[int, Array]] = None,
    dt: float = 1e-3,
    save_every: int = 1,
    seed: int = 1234,
    precomputed: Optional[Dict[str, Any]] = None,
) -> TrajectoryResult:
    """
    Quantum-trajectory simulation in the direct-sum strong-symmetry basis.

    Parameters
    ----------
    N
        Total number of atoms.
    gamma
        Collective decay rate Gamma in the paper.
    phases
        Piecewise-constant protocol stages. For your requested three-stage run,
        use e.g.
            [Phase(T1, Omega0, 0.0, 'phase1'),
             Phase(T2, Omega0, delta0, 'phase2'),
             Phase(T3, 0.0, 0.0, 'phase3')]
    sector_coeffs
        Initial amplitudes of the Nj sectors. Example:
            {N//2: 1.0, N//2 - 1: 1.0}
    internal_sector_states
        Optional internal |n_e>-basis vectors for individual sectors. If omitted,
        every populated sector starts in |n_e=0> (all active atoms in |down>).
    dt
        Base time step used for jump detection and propagation.
    save_every
        Save one snapshot every `save_every` accepted time steps.
    rng
        Optional NumPy random number generator.

    Returns
    -------
    TrajectoryResult
        Final wavefunction by sector, sector multiplicities, snapshots, and jump log.
    """
    if N <= 0:
        raise ValueError("N must be positive.")
    if gamma < 0.0:
        raise ValueError("gamma must be non-negative.")
    if dt <= 0.0:
        raise ValueError("dt must be positive.")
    if save_every <= 0:
        raise ValueError("save_every must be >= 1.")

    seed_seq = np.random.SeedSequence(seed).spawn(1)[0]
    rng = np.random.default_rng(seed_seq)

    # Keep dictionary inputs/outputs for compatibility, but use aligned lists in the hot loop.
    if precomputed is None:
        precomputed = build_precomputed_trajectory_data(N, gamma, phases, sector_coeffs, dt)

    sector_list = precomputed["sector_list"]
    ops_list = precomputed["ops_list"]
    multiplicities = precomputed["multiplicities"]
    dims = precomputed["dims"]
    phase_generators = precomputed["phase_generators"]
    phase_propagators = precomputed["phase_propagators"]

    initial_blocks = build_initial_sector_state(N, sector_coeffs, internal_sector_states)
    psi_blocks = renormalize_psi_blocks([initial_blocks[Nj] for Nj in sector_list])

    snapshots: List[TrajectorySnapshot] = [
        TrajectorySnapshot(
            time=0.0,
            sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
            norm=1.0,
            phase_index=0,
        )
    ]

    jump_times: List[float] = []
    current_time = 0.0
    accepted_steps = 0
    threshold = rng.random()

    for phase_index, phase in enumerate(phases):
        if phase.duration < 0.0:
            raise ValueError("Phase durations must be non-negative.")
        if phase.duration == 0.0:
            continue

        generators_list = phase_generators[phase_index]
        full_step_propagators = phase_propagators[phase_index]

        phase_end = current_time + phase.duration
        while current_time < phase_end - 1e-15:
            step = min(dt, phase_end - current_time)
            if abs(step - dt) <= 1e-15:
                trial = propagate_blocks_with_propagators(psi_blocks, full_step_propagators)
            else:
                trial = propagate_blocks(psi_blocks, generators_list, step)
            trial_norm2 = total_norm2_list(trial)

            if trial_norm2 > threshold:
                psi_blocks = trial
                current_time += step
                accepted_steps += 1

                if accepted_steps % save_every == 0:
                    snapshots.append(
                        TrajectorySnapshot(
                            time=current_time,
                            sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
                            norm=np.sqrt(trial_norm2),
                            phase_index=phase_index,
                        )
                    )
                continue

            lo, hi = 0.0, step
            pre_blocks = [psi.copy() for psi in psi_blocks]
            for _ in range(5):
                mid = 0.5 * (lo + hi)
                mid_state = propagate_blocks(pre_blocks, generators_list, mid)
                if total_norm2_list(mid_state) > threshold:
                    lo = mid
                else:
                    hi = mid

            tau = hi
            psi_blocks = propagate_blocks(pre_blocks, generators_list, tau)
            current_time += tau

            psi_blocks = renormalize_psi_blocks(psi_blocks)
            psi_blocks = apply_jump(psi_blocks, ops_list)
            jump_times.append(current_time)
            threshold = rng.random()

            snapshots.append(
                TrajectorySnapshot(
                    time=current_time,
                    sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
                    norm=1.0,
                    phase_index=phase_index,
                )
            )

            remainder = step - tau
            if remainder > 1e-15:
                trial = propagate_blocks(psi_blocks, generators_list, remainder)
                trial_norm2 = total_norm2_list(trial)
                if trial_norm2 <= threshold:
                    pass
                else:
                    psi_blocks = trial
                    current_time += remainder
                    accepted_steps += 1
                    if accepted_steps % save_every == 0:
                        snapshots.append(
                            TrajectorySnapshot(
                                time=current_time,
                                sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
                                norm=np.sqrt(trial_norm2),
                                phase_index=phase_index,
                            )
                        )

    psi_blocks = renormalize_psi_blocks(psi_blocks)
    if not snapshots or abs(snapshots[-1].time - current_time) > 1e-12:
        snapshots.append(
            TrajectorySnapshot(
                time=current_time,
                sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
                norm=1.0,
                phase_index=max(len(phases) - 1, 0),
            )
        )

    return TrajectoryResult(
        N=N,
        gamma=gamma,
        sectors=sector_list,
        sector_multiplicities=multiplicities,
        final_sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
        snapshots=snapshots,
        jump_times=jump_times,
        jump_count=len(jump_times),
        sector_dimensions=dims,
    )
