from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.j_moments import JMomentSeries
from parser.s_moments import SMomentSeries


def _validated_group_angles(
    moments: JMomentSeries | SMomentSeries,
    *,
    label: str,
    num_times: int,
) -> tuple[Array, Array]:
    """Return finite group-major angle arrays on the expected time grid."""
    if moments.theta is None or moments.phi is None:
        raise ValueError(f"{label} moments must contain group-resolved angles.")

    theta = np.asarray(moments.theta, dtype=float)
    phi = np.asarray(moments.phi, dtype=float)
    if theta.ndim != 2 or theta.shape != phi.shape:
        raise ValueError(
            f"{label} theta and phi must have matching shape "
            "(num_groups, num_times)."
        )
    if theta.shape[0] == 0 or theta.shape[1] != num_times:
        raise ValueError(
            f"{label} angles must contain at least one group and match the "
            "shared time grid."
        )
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(phi)):
        raise ValueError(f"{label} angles must be finite.")
    return theta, phi


def build_dressed_basis_transformation(
    j_moments: JMomentSeries,
    s_moments: SMomentSeries,
    *,
    tol: float = 1e-12,
) -> Array:
    """Build the time- and group-resolved ``(|c>, |j>, |s>)`` basis.

    The returned array has shape ``(num_times, num_groups, 3, 3)``. Its
    columns are the dressed states expressed in the lab basis
    ``(|u>, |d>, |e>)``. Thus ``U`` maps dressed coordinates to lab
    coordinates, while ``U.conj().swapaxes(-1, -2)`` maps lab coordinates to
    dressed ones for each time and group.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")

    t_j = np.asarray(j_moments.t, dtype=float)
    t_s = np.asarray(s_moments.t, dtype=float)
    if t_j.ndim != 1 or t_j.size == 0:
        raise ValueError("J moments must use a nonempty one-dimensional time grid.")
    if not np.all(np.isfinite(t_j)) or not np.all(np.isfinite(t_s)):
        raise ValueError("J and S time grids must be finite.")
    if t_s.shape != t_j.shape or not np.allclose(
        t_s,
        t_j,
        atol=tol,
        rtol=0.0,
    ):
        raise ValueError("J and S moments must share the same time grid.")

    theta_j, phi_j = _validated_group_angles(
        j_moments,
        label="J",
        num_times=t_j.size,
    )
    theta_s, phi_s = _validated_group_angles(
        s_moments,
        label="S",
        num_times=t_j.size,
    )
    if theta_s.shape != theta_j.shape:
        raise ValueError("J and S moments must contain the same physical groups.")

    angle_shape = theta_j.shape
    ket_one = np.stack(
        (
            np.zeros(angle_shape, dtype=np.complex128),
            np.cos(theta_j / 2.0),
            np.exp(-1j * phi_j) * np.sin(theta_j / 2.0),
        ),
        axis=-1,
    )
    ket_j = np.stack(
        (
            np.zeros(angle_shape, dtype=np.complex128),
            -np.sin(theta_j / 2.0),
            np.exp(-1j * phi_j) * np.cos(theta_j / 2.0),
        ),
        axis=-1,
    )

    ket_u = np.zeros(angle_shape + (3,), dtype=np.complex128)
    ket_u[..., 0] = 1.0
    phase_s = np.exp(-1j * phi_s)[..., None]
    ket_c = (
        np.cos(theta_s / 2.0)[..., None] * ket_u
        + phase_s * np.sin(theta_s / 2.0)[..., None] * ket_one
    )
    ket_s = (
        -np.sin(theta_s / 2.0)[..., None] * ket_u
        + phase_s * np.cos(theta_s / 2.0)[..., None] * ket_one
    )

    # Stack columns in the dressed coordinate order, then make time the first
    # public axis.
    transformation = np.swapaxes(
        np.stack((ket_c, ket_j, ket_s), axis=-1),
        0,
        1,
    ).copy()
    identity = np.eye(3, dtype=np.complex128)
    gram = np.einsum(
        "...ji,...jk->...ik",
        transformation.conjugate(),
        transformation,
        optimize=True,
    )
    if not np.allclose(gram, identity, atol=tol, rtol=0.0):
        raise ValueError("The constructed dressed-basis transformation is not unitary.")
    return transformation


__all__ = ["build_dressed_basis_transformation"]
