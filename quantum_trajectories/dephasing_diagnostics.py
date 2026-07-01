from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import numpy as np

from parser.common import ObservableSeries
from common.utils import active_manifold_angles, phase_change_times
from quantum_trajectories.aggregator import ensemble_observables, trajectory_observables
from parser.common import Array
from parser.quantum_trajectories import SectorKey, TrajectoryEnsemble, TrajectoryResult
from quantum_trajectories.state_helpers import total_norm2
from quantum_trajectories.utils import map_with_optional_pool


def _reference_result(result: Union[TrajectoryResult, TrajectoryEnsemble]) -> TrajectoryResult:
    if isinstance(result, TrajectoryEnsemble):
        if not result.trajectories:
            raise ValueError("TrajectoryEnsemble is empty.")
        return result.trajectories[0]
    return result


def _is_inhomogeneous_result(result: Union[TrajectoryResult, TrajectoryEnsemble]) -> bool:
    reference = _reference_result(result)
    return any(isinstance(sector_key, tuple) for sector_key in reference.sectors)


def _observable_series_for_result(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    tol: float,
    n_processes: Optional[int],
    averaged_observables: Optional[ObservableSeries],
) -> ObservableSeries:
    """
    Return trajectory or ensemble observables on the saved t_eval grid.

    For ensembles, this averages ordinary observables before any angle or
    effective-S construction, matching docs/instructions/bloch_vector_averaging.md.
    """
    if isinstance(result, TrajectoryEnsemble):
        reference = _reference_result(result)
        t_ref = np.asarray([snap.time for snap in reference.snapshots], dtype=float)
        if averaged_observables is not None:
            if len(averaged_observables.t) != len(t_ref) or not np.allclose(
                averaged_observables.t,
                t_ref,
                atol=1e-12,
                rtol=0.0,
            ):
                raise ValueError(
                    "averaged_observables must be defined on the same t_eval grid "
                    "as the ensemble trajectories."
                )
            return averaged_observables
        return ensemble_observables(result, tol=tol, n_processes=n_processes)

    if averaged_observables is not None:
        raise ValueError("averaged_observables is only valid for TrajectoryEnsemble input.")
    return trajectory_observables(result, tol=tol)


def _split_key(sector_key: SectorKey) -> tuple[int, ...]:
    if isinstance(sector_key, tuple):
        return tuple(int(v) for v in sector_key)
    return (int(sector_key),)


def _join_key(groups: Sequence[int]) -> SectorKey:
    if len(groups) == 1:
        return int(groups[0])
    return tuple(int(v) for v in groups)


def _block_inner_product(left: Mapping[SectorKey, Array], right: Mapping[SectorKey, Array]) -> complex:
    """
    Inner product between two direct-sum states with matching sector keys.
    """
    shared = set(left).intersection(right)
    return sum(np.vdot(left[key], right[key]) for key in shared)


def _add_to_block(out: dict[SectorKey, Array], sector_key: SectorKey, values: Array) -> None:
    if values.size == 0:
        return
    if sector_key in out:
        out[sector_key] = out[sector_key] + values
    else:
        out[sector_key] = np.asarray(values, dtype=np.complex128).copy()


def _apply_group_transition(
    blocks: Mapping[SectorKey, Array],
    *,
    group_sizes: Sequence[int],
    group_index: int,
    mu: str,
    nu: str,
) -> dict[SectorKey, Array]:
    """
    Apply A_{mu,nu}^{(group)} = sum_i |mu_i><nu_i| in one physical group.

    The same routine handles homogeneous keys (N_J) and inhomogeneous tuple
    keys (N_{J,1}, N_{J,2}).  It intentionally works at post-processing time
    only; no MCWF propagation data structures are changed.
    """
    out: dict[SectorKey, Array] = {}

    for sector_key, psi_flat in blocks.items():
        groups = list(_split_key(sector_key))
        if len(groups) != len(group_sizes):
            raise ValueError(
                f"Sector key {sector_key!r} is incompatible with group_sizes={group_sizes}."
            )
        if not (0 <= group_index < len(groups)):
            raise ValueError(f"group_index={group_index} is out of range for {len(groups)} groups.")

        dims = tuple(nj + 1 for nj in groups)
        psi = np.asarray(psi_flat, dtype=np.complex128).reshape(dims)
        target_blocks: dict[SectorKey, Array] = {}

        for idx in np.ndindex(*dims):
            amp = psi[idx]
            if amp == 0:
                continue

            groups_target = list(groups)
            idx_target = list(idx)
            nj_g = groups[group_index]
            n_e = idx[group_index]
            n_u = int(group_sizes[group_index]) - nj_g
            n_d = nj_g - n_e
            coeff = 0.0

            if mu == "u" and nu == "u":
                coeff = float(n_u)
            elif mu == "d" and nu == "d":
                coeff = float(n_d)
            elif mu == "e" and nu == "e":
                coeff = float(n_e)
            elif mu == "d" and nu == "e":
                if n_e <= 0:
                    continue
                coeff = np.sqrt(float(n_e * (n_d + 1)))
                idx_target[group_index] = n_e - 1
            elif mu == "e" and nu == "d":
                if n_d <= 0:
                    continue
                coeff = np.sqrt(float(n_d * (n_e + 1)))
                idx_target[group_index] = n_e + 1
            elif mu == "u" and nu == "d":
                if n_d <= 0:
                    continue
                coeff = np.sqrt(float((n_u + 1) * n_d))
                groups_target[group_index] = nj_g - 1
            elif mu == "d" and nu == "u":
                if n_u <= 0:
                    continue
                coeff = np.sqrt(float(n_u * (n_d + 1)))
                groups_target[group_index] = nj_g + 1
            elif mu == "u" and nu == "e":
                if n_e <= 0:
                    continue
                coeff = np.sqrt(float((n_u + 1) * n_e))
                groups_target[group_index] = nj_g - 1
                idx_target[group_index] = n_e - 1
            elif mu == "e" and nu == "u":
                if n_u <= 0:
                    continue
                coeff = np.sqrt(float(n_u * (n_e + 1)))
                groups_target[group_index] = nj_g + 1
                idx_target[group_index] = n_e + 1
            else:
                raise ValueError(f"Unsupported transition A_{{{mu}{nu}}}.")

            target_key = _join_key(groups_target)
            target_dims = tuple(nj + 1 for nj in groups_target)
            if target_key not in target_blocks:
                target_blocks[target_key] = np.zeros(target_dims, dtype=np.complex128)
            target_blocks[target_key][tuple(idx_target)] += coeff * amp

        for target_key, values in target_blocks.items():
            _add_to_block(out, target_key, values.reshape(-1))

    return out


def _apply_group_one_body_operator(
    blocks: Mapping[SectorKey, Array],
    *,
    group_sizes: Sequence[int],
    group_index: int,
    operator: Array,
) -> dict[SectorKey, Array]:
    """
    Apply one collective one-body operator in a selected physical group.
    """
    operator = np.asarray(operator, dtype=np.complex128)
    if operator.shape != (3, 3):
        raise ValueError(f"Expected a 3x3 single-particle operator, got {operator.shape}.")

    labels = ("u", "d", "e")
    out: dict[SectorKey, Array] = {}
    for row, mu in enumerate(labels):
        for col, nu in enumerate(labels):
            coeff = operator[row, col]
            if abs(coeff) <= 1e-15:
                continue
            transitioned = _apply_group_transition(
                blocks,
                group_sizes=group_sizes,
                group_index=group_index,
                mu=mu,
                nu=nu,
            )
            for target_key, values in transitioned.items():
                _add_to_block(out, target_key, coeff * values)
    return out


def _expect_group_one_body_operator(
    blocks: Mapping[SectorKey, Array],
    *,
    group_sizes: Sequence[int],
    group_index: int,
    operator: Array,
    norm2: float,
) -> complex:
    if norm2 <= 1e-15:
        return 0.0j
    applied = _apply_group_one_body_operator(
        blocks,
        group_sizes=group_sizes,
        group_index=group_index,
        operator=operator,
    )
    return _block_inner_product(blocks, applied) / norm2


def _dressed_one(theta: float, phi: float) -> Array:
    return np.array(
        [0.0, np.cos(theta / 2.0), np.exp(-1j * phi) * np.sin(theta / 2.0)],
        dtype=np.complex128,
    )


def _s_spin_operators(theta: float, phi: float) -> tuple[Array, Array, Array]:
    """
    Return Sx, Sy, Sz for the effective |u>, |1(theta, phi)> pseudospin.
    """
    ket_u = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
    ket_one = _dressed_one(theta, phi)
    sx_op = 0.5 * (
        np.outer(ket_one, ket_u.conjugate())
        + np.outer(ket_u, ket_one.conjugate())
    )
    sy_op = (
        np.outer(ket_one, ket_u.conjugate())
        - np.outer(ket_u, ket_one.conjugate())
    ) / (2.0j)
    sz_op = 0.5 * (
        np.outer(ket_u, ket_u.conjugate())
        - np.outer(ket_one, ket_one.conjugate())
    )
    return sx_op, sy_op, sz_op


def _group_angles_from_observables(
    obs: ObservableSeries,
    *,
    inhomogeneous: bool,
    tol: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build ensemble-level theta/phi series for each effective-S group.
    """
    if not inhomogeneous:
        return (
            np.asarray([obs.theta], dtype=float),
            np.asarray([obs.phi], dtype=float),
        )

    if (
        obs.Jx_groups is None
        or obs.Jy_groups is None
        or obs.Jz_groups is None
        or obs.N_e_groups is None
    ):
        raise ValueError("Inhomogeneous dephasing diagnostics require group-resolved observables.")

    theta_groups = []
    phi_groups = []
    for jx, jy, jz, ne in zip(obs.Jx_groups, obs.Jy_groups, obs.Jz_groups, obs.N_e_groups):
        theta, phi, _, _, _, _ = active_manifold_angles(jx, jy, jz, ne, tol=tol)
        theta_groups.append(theta)
        phi_groups.append(phi)

    return np.asarray(theta_groups, dtype=float), np.asarray(phi_groups, dtype=float)


def _group_sizes_for_result(result: Union[TrajectoryResult, TrajectoryEnsemble]) -> tuple[int, ...]:
    reference = _reference_result(result)
    if _is_inhomogeneous_result(result):
        if reference.N1 is None or reference.N2 is None:
            raise ValueError("Inhomogeneous results must store N1 and N2.")
        return int(reference.N1), int(reference.N2)
    return (int(reference.N),)


def _trajectory_s_components(
    traj: TrajectoryResult,
    *,
    t_ref: Array,
    theta_groups: Array,
    phi_groups: Array,
    group_sizes: Sequence[int],
    tol: float,
) -> np.ndarray:
    """
    Return raw <Sx>, <Sy>, <Sz> components for each group and saved time.
    """
    t = np.asarray([snap.time for snap in traj.snapshots], dtype=float)
    if len(t) != len(t_ref) or not np.allclose(t, t_ref, atol=1e-12, rtol=0.0):
        raise ValueError(
            "All trajectories must share the same t_eval snapshot grid. "
            "Run the ensemble through the common num_snapshots API."
        )

    group_count = len(group_sizes)
    components = np.zeros((group_count, len(t_ref), 3), dtype=float)

    for k, snap in enumerate(traj.snapshots):
        norm2 = total_norm2(snap.sector_blocks)
        if norm2 <= tol:
            continue

        for group_index in range(group_count):
            sx_op, sy_op, sz_op = _s_spin_operators(
                float(theta_groups[group_index, k]),
                float(phi_groups[group_index, k]),
            )
            components[group_index, k, 0] = float(
                np.real(
                    _expect_group_one_body_operator(
                        snap.sector_blocks,
                        group_sizes=group_sizes,
                        group_index=group_index,
                        operator=sx_op,
                        norm2=norm2,
                    )
                )
            )
            components[group_index, k, 1] = float(
                np.real(
                    _expect_group_one_body_operator(
                        snap.sector_blocks,
                        group_sizes=group_sizes,
                        group_index=group_index,
                        operator=sy_op,
                        norm2=norm2,
                    )
                )
            )
            components[group_index, k, 2] = float(
                np.real(
                    _expect_group_one_body_operator(
                        snap.sector_blocks,
                        group_sizes=group_sizes,
                        group_index=group_index,
                        operator=sz_op,
                        norm2=norm2,
                    )
                )
            )

    return components


def _trajectory_s_components_worker(args: tuple) -> np.ndarray:
    traj, t_ref, theta_groups, phi_groups, group_sizes, tol = args
    return _trajectory_s_components(
        traj,
        t_ref=t_ref,
        theta_groups=theta_groups,
        phi_groups=phi_groups,
        group_sizes=group_sizes,
        tol=tol,
    )


def dephasing_bloch_lengths(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    normalize: bool = True,
    tol: float = 1e-12,
    n_processes: Optional[int] = None,
    averaged_observables: Optional[ObservableSeries] = None,
) -> dict:
    """
    Compute effective-S Bloch-vector lengths for dephasing diagnostics.

    For ensembles, this averages <Sx>, <Sy>, and <Sz> over trajectories first,
    then computes the vector length. It deliberately does not average
    per-trajectory lengths.
    """
    reference = _reference_result(result)
    inhomogeneous = _is_inhomogeneous_result(result)
    group_sizes = _group_sizes_for_result(result)
    t = np.asarray([snap.time for snap in reference.snapshots], dtype=float)

    obs = _observable_series_for_result(
        result,
        tol=tol,
        n_processes=n_processes,
        averaged_observables=averaged_observables,
    )
    theta_groups, phi_groups = _group_angles_from_observables(
        obs,
        inhomogeneous=inhomogeneous,
        tol=tol,
    )

    if isinstance(result, TrajectoryEnsemble):
        per_traj_components = map_with_optional_pool(
            _trajectory_s_components_worker,
            [
                (traj, t, theta_groups, phi_groups, group_sizes, tol)
                for traj in result.trajectories
            ],
            n_processes=n_processes,
            progress_desc="dephasing S-Bloch components",
        )
        group_components = np.mean(np.asarray(per_traj_components, dtype=float), axis=0)
    else:
        group_components = _trajectory_s_components(
            result,
            t_ref=t,
            theta_groups=theta_groups,
            phi_groups=phi_groups,
            group_sizes=group_sizes,
            tol=tol,
        )

    total_components = np.sum(group_components, axis=0)
    group_lengths = np.linalg.norm(group_components, axis=2)
    total_length = np.linalg.norm(total_components, axis=1)
    total_normalization = float(reference.N) / 2.0
    group_normalizations = np.asarray(group_sizes, dtype=float) / 2.0

    if normalize:
        # Use the fixed physical maximum spin length, not the initial length
        # or the instantaneous active population. This makes dephasing curves
        # comparable across homogeneous and group-resolved inhomogeneous runs.
        total_length_plot = (
            total_length / total_normalization
            if total_normalization > tol
            else np.full_like(total_length, np.nan)
        )
        group_lengths_plot = np.full_like(group_lengths, np.nan)
        for group_index, scale in enumerate(group_normalizations):
            if scale > tol:
                group_lengths_plot[group_index] = group_lengths[group_index] / scale
        ylabel = r"$|\langle\mathbf{S}\rangle|/S_{\max}$"
    else:
        total_length_plot = total_length
        group_lengths_plot = group_lengths
        ylabel = r"$|\langle\mathbf{S}\rangle|$"

    delta_phi = None
    if inhomogeneous and theta_groups.shape[0] == 2:
        delta_phi = np.unwrap(phi_groups[0]) - np.unwrap(phi_groups[1])

    return {
        "t": t,
        "S_components_total": total_components,
        "S_components_groups": group_components,
        "length_total": total_length,
        "length_groups": group_lengths,
        "length_total_plot": total_length_plot,
        "length_groups_plot": group_lengths_plot,
        "theta_groups": theta_groups,
        "phi_groups": phi_groups,
        "delta_phi": delta_phi,
        "group_sizes": tuple(int(v) for v in group_sizes),
        "normalization_total": total_normalization if normalize else 1.0,
        "normalization_groups": (
            tuple(float(v) for v in group_normalizations)
            if normalize
            else tuple(1.0 for _ in group_sizes)
        ),
        "normalized": bool(normalize),
        "ylabel": ylabel,
        "inhomogeneous": bool(inhomogeneous),
    }


def plot_dephasing_bloch_lengths(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    normalize: bool = True,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    tol: float = 1e-12,
    n_processes: Optional[int] = None,
    averaged_observables: Optional[ObservableSeries] = None,
    plot_relative_phase: bool = True,
):
    """
    Plot effective-S Bloch-vector lengths as a standalone dephasing diagnostic.
    """
    import matplotlib.pyplot as plt

    data = dephasing_bloch_lengths(
        result,
        normalize=normalize,
        tol=tol,
        n_processes=n_processes,
        averaged_observables=averaged_observables,
    )
    reference = _reference_result(result)
    show_phase_panel = bool(plot_relative_phase and data["inhomogeneous"] and data["delta_phi"] is not None)

    if axes is None:
        if show_phase_panel:
            fig, axes_arr = plt.subplots(1, 2, figsize=(12, 4), sharex=True)
        else:
            fig, axes_arr = plt.subplots(1, 1, figsize=(7, 4))
            axes_arr = np.asarray([axes_arr])
    else:
        axes_arr = np.asarray(axes).ravel()
        expected = 2 if show_phase_panel else 1
        if axes_arr.size != expected:
            raise ValueError(f"Expected {expected} axes, got {axes_arr.size}.")
        fig = axes_arr[0].figure

    t = data["t"]
    length_ax = axes_arr[0]
    length_ax.plot(t, data["length_total_plot"], linewidth=1.8, label="total")
    if data["inhomogeneous"]:
        for group_index, values in enumerate(data["length_groups_plot"], start=1):
            length_ax.plot(t, values, linewidth=1.5, label=f"group {group_index}")

    if len(reference.phases) >= 2:
        for phase_time in phase_change_times(reference.phases):
            for ax in axes_arr:
                ax.axvline(phase_time, linestyle="--", color="black", alpha=0.45)

    length_ax.set_xlabel(r"$\Gamma t$")
    length_ax.set_ylabel(data["ylabel"])
    length_ax.set_title("Dephasing diagnostic")
    length_ax.grid(alpha=0.3)
    length_ax.legend()
    length_ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)

    if show_phase_panel:
        phase_ax = axes_arr[1]
        phase_ax.plot(t, data["delta_phi"], linewidth=1.8, label=r"$\Delta\phi=\phi_1-\phi_2$")
        phase_ax.set_xlabel(r"$\Gamma t$")
        phase_ax.set_ylabel(r"$\Delta\phi$")
        phase_ax.set_title("Relative group phase")
        phase_ax.grid(alpha=0.3)
        phase_ax.legend()
        phase_ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)

    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return data, fig, axes_arr


def plot_compare_dephasing_bloch_lengths(
    homogeneous_result: Union[TrajectoryResult, TrajectoryEnsemble],
    inhomogeneous_result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    normalize: bool = True,
    ax=None,
    output_path: Optional[Union[str, Path]] = None,
    tol: float = 1e-12,
    n_processes: Optional[int] = None,
    homogeneous_averaged_observables: Optional[ObservableSeries] = None,
    inhomogeneous_averaged_observables: Optional[ObservableSeries] = None,
):
    """
    Compute and plot homogeneous-vs-inhomogeneous dephasing diagnostics.

    This mirrors the comparison notebook use case: the homogeneous total length,
    the inhomogeneous total length, and the inhomogeneous group-resolved lengths
    are shown on one axis.  The dephasing lengths are computed inside this
    helper so callers do not need to duplicate plotting logic.
    """
    import matplotlib.pyplot as plt

    homogeneous_data = dephasing_bloch_lengths(
        homogeneous_result,
        normalize=normalize,
        tol=tol,
        n_processes=n_processes,
        averaged_observables=homogeneous_averaged_observables,
    )
    inhomogeneous_data = dephasing_bloch_lengths(
        inhomogeneous_result,
        normalize=normalize,
        tol=tol,
        n_processes=n_processes,
        averaged_observables=inhomogeneous_averaged_observables,
    )

    if ax is None:
        fig, axis = plt.subplots(figsize=(8, 4))
    else:
        axis = ax
        fig = axis.figure

    axis.plot(
        homogeneous_data["t"],
        homogeneous_data["length_total_plot"],
        linewidth=2.0,
        label="homogeneous total",
    )
    axis.plot(
        inhomogeneous_data["t"],
        inhomogeneous_data["length_total_plot"],
        linewidth=2.0,
        label="inhomogeneous total",
    )
    if inhomogeneous_data["inhomogeneous"]:
        for group_index, group_length in enumerate(
            inhomogeneous_data["length_groups_plot"],
            start=1,
        ):
            axis.plot(
                inhomogeneous_data["t"],
                group_length,
                linestyle="--",
                linewidth=1.4,
                alpha=0.8,
                label=f"inhomogeneous group {group_index}",
            )

    reference = _reference_result(inhomogeneous_result)
    if len(reference.phases) >= 2:
        for phase_time in phase_change_times(reference.phases):
            axis.axvline(phase_time, linestyle="--", color="black", alpha=0.35)

    axis.set_xlabel(r"$\Gamma t$")
    axis.set_ylabel(homogeneous_data["ylabel"])
    axis.set_title("Dephasing Bloch-vector length")
    axis.grid(alpha=0.3)
    axis.legend()
    axis.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return {
        "homogeneous": homogeneous_data,
        "inhomogeneous": inhomogeneous_data,
    }, fig, axis
