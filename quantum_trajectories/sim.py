from __future__ import annotations

from quantum_trajectories.parser import (
    Array,
    Phase,
    SectorOperators,
    TrajectorySnapshot,
    TrajectoryResult,
)
from quantum_trajectories.operator_helpers import (
    build_sector_ops,
    sector_multiplicity,
)
from quantum_trajectories.state_helpers import (
    build_initial_sector_state,
    total_norm2,
)

from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply


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



def renormalize_blocks(blocks: Dict[int, Array]) -> Dict[int, Array]:
    nrm = np.sqrt(total_norm2(blocks))
    if nrm == 0.0:
        raise RuntimeError("Wavefunction has zero norm.")
    return {Nj: psi / nrm for Nj, psi in blocks.items()}


def propagate_blocks(
    blocks: Mapping[int, Array],
    generators: Mapping[int, csc_matrix],
    dt: float,
) -> Dict[int, Array]:
    """
    Propagate each sector block forward in time by dt under its effective generator using expm_multiply.
    """
    out: Dict[int, Array] = {}
    for Nj, psi in blocks.items():
        if psi.size == 0:
            out[Nj] = psi.copy()
            continue
        out[Nj] = expm_multiply((-1j * generators[Nj]) * dt, psi)
    return out



def apply_jump(blocks: Mapping[int, Array], ops_map: Mapping[int, SectorOperators]) -> Dict[int, Array]:
    jumped = {Nj: jump_for_sector(ops_map[Nj], psi) for Nj, psi in blocks.items()}
    if total_norm2(jumped) == 0.0:
        # No further jumps possible; keep the pre-jump state instead of crashing.
        return {Nj: psi.copy() for Nj, psi in blocks.items()}
    return renormalize_blocks(jumped)


# -----------------------------------------------------------------------------
# Main Monte Carlo trajectory engine
# -----------------------------------------------------------------------------


def simulate_nj_sector_trajectory(
    N: int,
    gamma: float,
    phases: Sequence[Phase],
    sector_coeffs: Mapping[int, complex],
    *,
    internal_sector_states: Optional[Mapping[int, Array]] = None,
    dt: float = 1e-3,
    save_every: int = 1,
    seed: int = 1234,
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

    rng=np.random.default_rng(seed)

    # Sorted dictionaty of key: sector Nj, value: coefficient for that sector
    sectors = sorted(sector_coeffs.keys())
    # Dictionary of key: sector Nj, value: operators (Jp,Jm...) for that sector
    ops_map = {Nj: build_sector_ops(Nj) for Nj in sectors}
    # Dictionary of key: sector Nj, value: (number of states with that total excitation number / Nj
    multiplicities = {Nj: sector_multiplicity(N, Nj) for Nj in sectors}
    # Dictionary of key: sector Nj, value: dimension of that sector's Hilbert space in the e,d manifold (Nj + 1)
    dims = {Nj: Nj + 1 for Nj in sectors}

    # Dictionary of key: sector Nj, value: normalized state vector in that sector's symmetric |n_e> basis (1,0,0,...,0), multiplied by the normalized sector coefficient.
    # blocks = {Nj: (1,0,0,...,0) * coeff for Nj, coeff in coeffs.items()}
    # (1,0,0,...,0) is the down state in each sector
    blocks = build_initial_sector_state(N, sector_coeffs, internal_sector_states)
    blocks = renormalize_blocks(blocks)

    # Create initial snapshot at time 0 with the initial state.
    snapshots: List[TrajectorySnapshot] = [
        TrajectorySnapshot(
            time=0.0,
            sector_blocks={Nj: psi.copy() for Nj, psi in blocks.items()},
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

        # Dictionary of key: sector Nj, value: effective non-Hermitian Hamiltonian for that sector in this phase.
        generators = {
            Nj: heff_for_sector(ops_map[Nj], phase.omega, phase.delta, gamma)
            for Nj in sectors
        }

        phase_end = current_time + phase.duration
        while current_time < phase_end - 1e-15:
            step = min(dt, phase_end - current_time)
            # Propagate state forward by non-unitary effective Hamiltonian
            trial = propagate_blocks(blocks, generators, step)
            trial_norm2 = total_norm2(trial)

            # No jump
            if trial_norm2 > threshold:
                blocks = trial
                current_time += step
                accepted_steps += 1

                if accepted_steps % save_every == 0:
                    snapshots.append(
                        TrajectorySnapshot(
                            time=current_time,
                            sector_blocks={Nj: psi.copy() for Nj, psi in blocks.items()},
                            norm=np.sqrt(trial_norm2),
                            phase_index=phase_index,
                        )
                    )
                continue

            # A jump happened inside this step. Refine the jump time by bisection.
            lo, hi = 0.0, step
            pre = {Nj: psi.copy() for Nj, psi in blocks.items()}
            for _ in range(15):
                mid = 0.5 * (lo + hi)
                mid_state = propagate_blocks(pre, generators, mid)
                if total_norm2(mid_state) > threshold:
                    lo = mid
                else:
                    hi = mid

            tau = hi
            blocks = propagate_blocks(pre, generators, tau)
            current_time += tau

            # Normalize pre-jump state, apply the regular jump operator, renormalize.
            blocks = renormalize_blocks(blocks)
            blocks = apply_jump(blocks, ops_map)
            jump_times.append(current_time)
            threshold = rng.random()

            snapshots.append(
                TrajectorySnapshot(
                    time=current_time,
                    sector_blocks={Nj: psi.copy() for Nj, psi in blocks.items()},
                    norm=1.0,
                    phase_index=phase_index,
                )
            )

            # Continue with the remainder of the step in the same phase.
            remainder = step - tau
            if remainder > 1e-15:
                trial = propagate_blocks(blocks, generators, remainder)
                trial_norm2 = total_norm2(trial)
                if trial_norm2 <= threshold:
                    # If another jump happens immediately, let the next loop catch it.
                    # Using a smaller dt is the intended way to resolve such clustering.
                    pass
                else:
                    blocks = trial
                    current_time += remainder
                    accepted_steps += 1
                    if accepted_steps % save_every == 0:
                        snapshots.append(
                            TrajectorySnapshot(
                                time=current_time,
                                sector_blocks={Nj: psi.copy() for Nj, psi in blocks.items()},
                                norm=np.sqrt(trial_norm2),
                                phase_index=phase_index,
                            )
                        )

    # Return the physical state normalized to 1.
    blocks = renormalize_blocks(blocks)
    if not snapshots or abs(snapshots[-1].time - current_time) > 1e-12:
        snapshots.append(
            TrajectorySnapshot(
                time=current_time,
                sector_blocks={Nj: psi.copy() for Nj, psi in blocks.items()},
                norm=1.0,
                phase_index=max(len(phases) - 1, 0),
            )
        )

    return TrajectoryResult(
        N=N,
        gamma=gamma,
        sectors=sectors,
        sector_multiplicities=multiplicities,
        final_sector_blocks={Nj: psi.copy() for Nj, psi in blocks.items()},
        snapshots=snapshots,
        jump_times=jump_times,
        jump_count=len(jump_times),
        sector_dimensions=dims,
    )
