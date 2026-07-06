from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
from scipy.sparse import csc_matrix, eye
from scipy.sparse.linalg import expm, expm_multiply

from solvers.mcwf.operator_helpers import (
    build_sector_ops_for_key,
    sector_multiplicity,
    split_sector_key,
    total_active_atoms_in_sector,
    two_group_sector_multiplicity,
)
from parser.common import Array, Phase
from common.utils.phases import phase_boundary_times
from parser.mcwf import (
    SectorKey,
    SectorOperators,
    TrajectoryResult,
    TrajectorySnapshot,
)
from solvers.mcwf.state_helpers import (
    build_initial_sector_state,
    total_norm2,
)


# -----------------------------------------------------------------------------
# Effective generators and propagation
# -----------------------------------------------------------------------------


def build_t_eval_from_phases(phases: Sequence[Phase], num_snapshots: int) -> np.ndarray:
    """
    Construct the common saved-time grid for MCWF trajectories.

    All trajectories in an ensemble use the same t_eval values so observables
    and squeezing moments can be averaged at identical physical times.
    """
    if num_snapshots < 2:
        raise ValueError("num_snapshots must be at least 2.")
    total_time = float(phase_boundary_times(phases)[-1])
    return np.linspace(0.0, total_time, num_snapshots, dtype=float)


def build_phase_jump_operator_for_sector(
    ops: SectorOperators,
    omega: float,
    Gamma: float,
    *,
    shifted_jump_operator: bool = False,
) -> csc_matrix:
    """
    Build the reduced-basis jump operator for one sector and one protocol phase.

    In the regular picture, the jump operator is simply
        l = J_-.

    In the shifted picture used for the cavity-output formulation, it becomes
        l = J_- + i Omega / Gamma.

    Here `ops` contains the reduced operators for a fixed Nj sector, so the
    returned matrix acts only within that sector's (Nj + 1)-dimensional basis.
    """
    unshifted_jump = ops.A_weighted if ops.A_weighted is not None else ops.Jm
    if not shifted_jump_operator:
        return unshifted_jump

    if Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    identity = eye(unshifted_jump.shape[0], format="csc", dtype=np.complex128)
    return (unshifted_jump + (1j * omega / Gamma) * identity).tocsc()


def heff_for_sector(
    ops: SectorOperators,
    omega: float,
    delta: float,
    Gamma: float,
    *,
    shifted_jump_operator: bool = False,
    jump_operator: Optional[csc_matrix] = None,
) -> csc_matrix:
    """
    Construct the effective Hamiltonian for a given sector.

    Regular H_delta and regular jump operator from the paper:
        H_delta = Omega J_x - delta N_e,
        l = J_-,
        H_eff = H_delta - i (Gamma/2) J_+ J_-.

    In MCQT, between jumps the state evolves with the non-Hermitian effective Hamiltonian
        H_eff = H_sys - i/2 sum_l l^dagger l
    For the homogeneous collective decay, the only jump operator is J_-, so the sum reduces to J_+ J_-.
        H_delta - i/2 Gamma J_+ J_-,
    """
    if not shifted_jump_operator:
        drive_op = ops.J_drive if ops.J_drive is not None else ops.J_x
        decay_term = ops.AdagA_weighted if ops.AdagA_weighted is not None else ops.JpJm
        H = omega * drive_op - delta * ops.N_e
        Heff = H - 0.5j * Gamma * decay_term
        return Heff.tocsc()

    if jump_operator is None:
        jump_operator = build_phase_jump_operator_for_sector(
            ops,
            omega,
            Gamma,
            shifted_jump_operator=True,
        )

    H = -delta * ops.N_e
    jump_dag_jump = (jump_operator.conjugate().transpose() @ jump_operator).tocsc()
    Heff = H - 0.5j * Gamma * jump_dag_jump
    return Heff.tocsc()


def jump_for_sector(jump_operator: csc_matrix, psi: Array) -> Array:
    return jump_operator.dot(psi)


def blocks_list_to_dict(sector_list: Sequence[SectorKey], psi_blocks: Sequence[Array]) -> Dict[SectorKey, Array]:
    return {Nj: psi.copy() for Nj, psi in zip(sector_list, psi_blocks)}


def total_norm2_list(psi_blocks: Sequence[Array]) -> float:
    return float(sum(np.vdot(psi, psi).real for psi in psi_blocks))


def renormalize_blocks(blocks: Dict[SectorKey, Array]) -> Dict[SectorKey, Array]:
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


def apply_jump(psi_blocks: Sequence[Array], jump_operators_list: Sequence[csc_matrix]) -> List[Array]:
    jumped = [jump_for_sector(jump_operator, psi) for jump_operator, psi in zip(jump_operators_list, psi_blocks)]
    if total_norm2_list(jumped) == 0.0:
        # No further jumps possible; keep the pre-jump state instead of crashing.
        return [psi.copy() for psi in psi_blocks]
    return renormalize_psi_blocks(jumped)


def build_precomputed_trajectory_data(
    Ni: Sequence[int],
    omega_i: Sequence[float],
    Gamma: float,
    phases: Sequence[Phase],
    sector_coeffs: Mapping[SectorKey, complex],
    dt: float,
    *,
    shifted_jump_operator: bool = False,
) -> Dict[str, Any]:
    if not Ni:
        raise ValueError("Ni must contain at least one group size.")
    if len(omega_i) != len(Ni):
        raise ValueError("omega_i must have the same length as Ni.")

    Ni = [int(group_size) for group_size in Ni]
    omega_i = [float(coupling) for coupling in omega_i]

    # Sorted list of populated strong-symmetry sectors, e.g. [N//2 - 1, N//2, N//2 + 1].
    sector_list = sorted(
        sector_coeffs.keys(),
        key=lambda key: split_sector_key(key),
    )

    # One reduced-basis operator bundle per sector:
    # [SectorOperators(Nj_0), SectorOperators(Nj_1), ...].
    ops_list = [
        build_sector_ops_for_key(
            sector_key,
            Ni=Ni,
            omega_i=omega_i,
        )
        for sector_key in sector_list
    ]
    # Sector multiplicity lookup, e.g. {500: ...} or {(250, 250): ...}.
    multiplicities = {
        sector_key: (
            sector_multiplicity(Ni[0], int(sector_key))
            if not isinstance(sector_key, tuple)
            else two_group_sector_multiplicity(Ni[0], Ni[1], int(sector_key[0]), int(sector_key[1]))
        )
        for sector_key in sector_list
    }
    # Reduced Hilbert-space dimension in each sector: Nj + 1 or (Nj1 + 1)(Nj2 + 1).
    dims = {sector_key: ops.Jm.shape[0] for sector_key, ops in zip(sector_list, ops_list)}
    # Phase- and sector-resolved jump operators:
    # phase_jump_operators[phase_index][sector_index] = l_{phase,Nj}.
    phase_jump_operators = [
        [
            build_phase_jump_operator_for_sector(
                ops,
                phase.omega,
                Gamma,
                shifted_jump_operator=shifted_jump_operator,
            )
            for ops in ops_list
        ]
        for phase in phases
    ]
    # Phase- and sector-resolved non-Hermitian generators:
    # phase_generators[phase_index][sector_index] = H_eff^{(phase,Nj)}.
    phase_generators = [
        [
            heff_for_sector(
                ops,
                phase.omega,
                phase.delta,
                Gamma,
                shifted_jump_operator=shifted_jump_operator,
                jump_operator=jump_operator,
            )
            for ops, jump_operator in zip(ops_list, jump_operators_list)
        ]
        for phase, jump_operators_list in zip(phases, phase_jump_operators)
    ]
    # Precomputed full-dt propagators:
    # phase_propagators[phase_index][sector_index] = exp(-i H_eff dt).
    phase_propagators = [
        [expm((-1j * generator) * dt).tocsc() for generator in generators_list]
        for generators_list in phase_generators
    ]

    return {
        "sector_list": sector_list,  # [Nj_0, Nj_1, ...]
        "ops_list": ops_list,  # [SectorOperators(Nj_0), SectorOperators(Nj_1), ...]
        "multiplicities": multiplicities,  # {Nj: sector multiplicity}
        "dims": dims,  # {Nj: reduced sector dimension = Nj + 1}
        "phase_jump_operators": phase_jump_operators,  # [phase][sector] -> l_{phase,Nj}
        "phase_generators": phase_generators,  # [phase][sector] -> H_eff^{(phase,Nj)}
        "phase_propagators": phase_propagators,  # [phase][sector] -> exp(-i H_eff dt)
    }


# -----------------------------------------------------------------------------
# Main Monte Carlo trajectory engine
# -----------------------------------------------------------------------------


def _simulate_single_trajectory(
    Ni: Sequence[int],
    omega_i: Sequence[float],
    Gamma: float,
    phases: Sequence[Phase],
    sector_coeffs: Mapping[SectorKey, complex],
    *,
    dt: float = 1e-3,
    t_eval: Array,
    seed_sequence: np.random.SeedSequence,
    shifted_jump_operator: bool = False,
    precomputed: Dict[str, Any],
) -> TrajectoryResult:
    """
    Quantum-trajectory simulation in the direct-sum strong-symmetry basis.

    Parameters
    ----------
    Ni
        Group sizes. Homogeneous runs should still pass one-entry lists.
    Gamma
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
    dt
        Base time step used for jump detection and propagation.
    t_eval
        Explicit saved-time grid shared by the ensemble. The simulator saves
        the state exactly at these times by splitting internal evolution steps
        when needed.
    omega_i
        Group couplings. Must have the same length as Ni.
    seed_sequence
        Child seed assigned by the ensemble runner for this one trajectory.

    Returns
    -------
    TrajectoryResult
        Final wavefunction by sector, sector multiplicities, snapshots, and jump log.
    """
    assert len(omega_i) == len(Ni)
    assert precomputed is not None

    N = sum(Ni)
    assert N > 0

    rng = np.random.default_rng(seed_sequence)

    sector_list = precomputed["sector_list"]
    multiplicities = precomputed["multiplicities"]
    dims = precomputed["dims"]
    phase_jump_operators = precomputed["phase_jump_operators"]
    phase_generators = precomputed["phase_generators"]
    phase_propagators = precomputed["phase_propagators"]

    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or t_eval.size < 2:
        raise ValueError("t_eval must be a one-dimensional array with at least two points.")
    if np.any(np.diff(t_eval) <= 0.0):
        raise ValueError("t_eval must be strictly increasing.")
    total_time = float(phase_boundary_times(phases)[-1])
    if abs(float(t_eval[0])) > 1e-12:
        raise ValueError("The first t_eval point must be 0.0.")
    if abs(float(t_eval[-1]) - total_time) > 1e-9:
        raise ValueError("The last t_eval point must match the total protocol time.")

    next_eval_idx = 1

    initial_blocks = build_initial_sector_state(N, sector_coeffs)
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
    threshold = rng.random()
    # Count actual propagation calls. This lets the runtime diagnostic include
    # bisection midpoint propagations and other variable-step fallbacks, not
    # only the outer-loop step attempts.
    total_step_count = 0
    non_precomputed_step_count = 0

    def maybe_save_snapshot() -> None:
        nonlocal next_eval_idx
        while next_eval_idx < len(t_eval) and current_time >= t_eval[next_eval_idx] - 1e-15:
            if abs(current_time - t_eval[next_eval_idx]) > 1e-9:
                raise RuntimeError(
                    "Internal MCWF evolution missed a requested t_eval output time. "
                    "This indicates the step-splitting logic failed."
                )
            snapshots.append(
                TrajectorySnapshot(
                    time=float(t_eval[next_eval_idx]),
                    sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
                    norm=np.sqrt(total_norm2_list(psi_blocks)),
                    phase_index=phase_index,
                )
            )
            next_eval_idx += 1

    for phase_index, phase in enumerate(phases):
        if phase.duration < 0.0:
            raise ValueError("Phase durations must be non-negative.")
        if phase.duration == 0.0:
            continue

        jump_operators_list = phase_jump_operators[phase_index]
        generators_list = phase_generators[phase_index]
        full_step_propagators = phase_propagators[phase_index]

        phase_end = current_time + phase.duration
        while current_time < phase_end - 1e-15:
            next_eval_time = t_eval[next_eval_idx] if next_eval_idx < len(t_eval) else np.inf
            step = min(dt, phase_end - current_time, next_eval_time - current_time)
            if abs(step - dt) <= 1e-15:
                total_step_count += 1
                trial = propagate_blocks_with_propagators(psi_blocks, full_step_propagators)
            else:
                # Any step that is shorter than the base dt cannot use the
                # precomputed full-step propagator and falls back to the
                # variable-step propagation path instead.
                total_step_count += 1
                non_precomputed_step_count += 1
                trial = propagate_blocks(psi_blocks, generators_list, step)
            trial_norm2 = total_norm2_list(trial)

            if trial_norm2 > threshold:
                psi_blocks = trial
                current_time += step
                maybe_save_snapshot()
                continue

            lo, hi = 0.0, step
            pre_blocks = [psi.copy() for psi in psi_blocks]
            for _ in range(5):
                mid = 0.5 * (lo + hi)
                total_step_count += 1
                non_precomputed_step_count += 1
                mid_state = propagate_blocks(pre_blocks, generators_list, mid)
                if total_norm2_list(mid_state) > threshold:
                    lo = mid
                else:
                    hi = mid

            tau = hi
            total_step_count += 1
            non_precomputed_step_count += 1
            psi_blocks = propagate_blocks(pre_blocks, generators_list, tau)
            current_time += tau

            psi_blocks = renormalize_psi_blocks(psi_blocks)
            psi_blocks = apply_jump(psi_blocks, jump_operators_list)
            jump_times.append(current_time)
            threshold = rng.random()
            maybe_save_snapshot()

            remainder = step - tau
            if remainder > 1e-15:
                total_step_count += 1
                non_precomputed_step_count += 1
                trial = propagate_blocks(psi_blocks, generators_list, remainder)
                trial_norm2 = total_norm2_list(trial)
                if trial_norm2 <= threshold:
                    pass
                else:
                    psi_blocks = trial
                    current_time += remainder
                    maybe_save_snapshot()

    psi_blocks = renormalize_psi_blocks(psi_blocks)
    if next_eval_idx < len(t_eval):
        if abs(current_time - t_eval[next_eval_idx]) > 1e-9:
            raise RuntimeError(
                "Final MCWF time does not match the last requested t_eval point."
            )
        snapshots.append(
            TrajectorySnapshot(
                time=float(t_eval[next_eval_idx]),
                sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
                norm=1.0,
                phase_index=max(len(phases) - 1, 0),
            )
        )
        next_eval_idx += 1

    if next_eval_idx != len(t_eval):
        raise RuntimeError("Did not save all requested MCWF t_eval snapshots.")

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
        final_sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
        snapshots=snapshots,
        jump_times=jump_times,
        jump_count=len(jump_times),
        total_step_count=total_step_count,
        non_precomputed_step_count=non_precomputed_step_count,
    )
