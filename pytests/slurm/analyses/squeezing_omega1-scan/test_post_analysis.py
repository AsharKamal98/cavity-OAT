from __future__ import annotations

import numpy as np

from parser.moments import MomentSeries
from post_analysis import (
    build_dressed_basis_transformation,
    compute_gell_mann_rotation,
    compute_j_moments,
    compute_mean_generator_moments,
    compute_s_moments,
    compute_squeezing,
    rotate_generator_moments,
)


def _compute_direct_moments(case) -> MomentSeries:
    moments = MomentSeries(
        t=case.ensemble.moments.t,
        metadata=case.simulation_metadata,
    )
    moments.generator_moments = case.ensemble.moments
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


def _assert_fields_allclose(actual, expected, field_names) -> None:
    for field_name in field_names:
        actual_value = getattr(actual, field_name)
        expected_value = getattr(expected, field_name)
        if isinstance(expected_value, tuple):
            assert len(actual_value) == len(expected_value)
            for actual_group, expected_group in zip(actual_value, expected_value):
                np.testing.assert_allclose(
                    actual_group,
                    expected_group,
                    atol=1e-12,
                    rtol=0.0,
                    equal_nan=True,
                )
        else:
            np.testing.assert_allclose(
                actual_value,
                expected_value,
                atol=1e-12,
                rtol=0.0,
                equal_nan=True,
            )


def _assert_moments_allclose(actual: MomentSeries, expected: MomentSeries) -> None:
    for field_name in ("generator_moments", "mean_generator_moments"):
        actual_generator_moments = getattr(actual, field_name)
        expected_generator_moments = getattr(expected, field_name)
        _assert_fields_allclose(
            actual_generator_moments,
            expected_generator_moments,
            ("t", "first_moments", "second_moments"),
        )

    _assert_fields_allclose(
        actual.J,
        expected.J,
        (
            "t",
            "integration_phase_index",
            "x",
            "y",
            "z",
            "length",
            "nx",
            "ny",
            "nz",
            "N_e",
            "N_j",
            "theta",
            "phi",
        ),
    )
    _assert_fields_allclose(
        actual.J.aggregated,
        expected.J.aggregated,
        (
            "x",
            "y",
            "z",
            "length",
            "nx",
            "ny",
            "nz",
            "N_e",
            "N_j",
            "theta",
            "phi",
        ),
    )
    _assert_fields_allclose(
        actual.S,
        expected.S,
        (
            "t",
            "integration_phase_index",
            "x",
            "y",
            "z",
            "length",
            "nx",
            "ny",
            "nz",
            "theta",
            "phi",
        ),
    )
    assert actual.S.basis is not None
    assert expected.S.basis is not None
    _assert_fields_allclose(
        actual.S.basis,
        expected.S.basis,
        ("zero_state", "one_state"),
    )

    assert len(actual.squeezing.group_resolved) == len(
        expected.squeezing.group_resolved
    )
    for actual_group, expected_group in zip(
        actual.squeezing.group_resolved,
        expected.squeezing.group_resolved,
    ):
        assert actual_group.N_atoms == expected_group.N_atoms
        _assert_fields_allclose(
            actual_group,
            expected_group,
            (
                "t",
                "covariance",
                "covariance_eigvals",
                "lambda_min",
                "xi2",
                "xi2_db",
                "N_c",
                "excited_fraction_active",
            ),
        )
    assert actual.squeezing.full_system is None
    assert expected.squeezing.full_system is None


def test_combined_post_analysis_matches_direct_full_ensemble(
    squeezing_post_analysis_case,
    squeezing_post_analysis_module,
):
    """
    Compare analysis after loading two chunks with direct full-ensemble analysis.

    Assert that raw, mean, J, S, and squeezing moments are unchanged.
    """
    case = squeezing_post_analysis_case
    combined, simulation_metadata = (
        squeezing_post_analysis_module._load_combined_ensemble(
            case.input_dir,
            case.filename,
            len(case.chunk_sizes),
        )
    )

    actual = squeezing_post_analysis_module._compute_moments(
        combined,
        simulation_metadata,
    )
    expected = _compute_direct_moments(case)

    _assert_moments_allclose(actual, expected)


def test_combined_ensemble_preserves_trajectory_order(
    squeezing_post_analysis_case,
    squeezing_post_analysis_module,
):
    """
    Load uneven trajectory chunks in array-index order.

    Assert that every moment row remains ordered and its seed records its chunk.
    """
    case = squeezing_post_analysis_case

    combined, _ = squeezing_post_analysis_module._load_combined_ensemble(
        case.input_dir,
        case.filename,
        len(case.chunk_sizes),
    )

    np.testing.assert_allclose(
        combined.moments.first_moments,
        case.ensemble.moments.first_moments,
        atol=0.0,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        combined.moments.second_moments,
        case.ensemble.moments.second_moments,
        atol=0.0,
        rtol=0.0,
    )
    assert combined.seeds == [(0, 0), (1, 1), (1, 2)]
