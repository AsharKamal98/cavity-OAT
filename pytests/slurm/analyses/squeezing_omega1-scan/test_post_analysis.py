from __future__ import annotations

import numpy as np

from parser.moments import MomentSeries
from post_analysis.squeezing import compute_squeezing
from solvers.mcwf.j_moments import compute_mcwf_j_moments
from solvers.mcwf.s_moments import compute_mcwf_s_moments
from solvers.mcwf.squeezing_moments import compute_mcwf_squeezing_moments


def _compute_direct_moments(case) -> MomentSeries:
    moments = MomentSeries(
        t=case.ensemble.metadata.t_eval,
        metadata=case.simulation_metadata,
    )
    moments.J = compute_mcwf_j_moments(case.ensemble, n_processes=1)
    moments.S = compute_mcwf_s_moments(
        case.ensemble,
        moments.J,
        n_processes=1,
    )
    moments.squeezing_moments = compute_mcwf_squeezing_moments(
        case.ensemble,
        moments.J,
        moments.S,
        n_processes=1,
    )
    moments.squeezing = compute_squeezing(moments.squeezing_moments)
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
            "jump_rate",
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

    assert len(actual.squeezing_moments) == len(expected.squeezing_moments)
    for actual_group, expected_group in zip(
        actual.squeezing_moments,
        expected.squeezing_moments,
    ):
        assert actual_group.N_group == expected_group.N_group
        _assert_fields_allclose(
            actual_group,
            expected_group,
            (
                "t",
                "mean_operators",
                "second_moments",
                "N_c",
                "excited_fraction_active",
                "theta_J",
                "phi_J",
                "theta_S",
                "phi_S",
            ),
        )

    assert len(actual.squeezing) == len(expected.squeezing)
    for actual_group, expected_group in zip(
        actual.squeezing,
        expected.squeezing,
    ):
        assert actual_group.N_group == expected_group.N_group
        _assert_fields_allclose(
            actual_group,
            expected_group,
            (
                "t",
                "xi2",
                "xi2_db",
                "lambda_min",
                "covariance_eigvals",
                "N_c",
                "excited_fraction_active",
                "theta_J",
                "phi_J",
                "theta_S",
                "phi_S",
            ),
        )


def test_combined_post_analysis_matches_direct_full_ensemble(
    squeezing_post_analysis_case,
    squeezing_post_analysis_module,
):
    """
    Compare analysis after loading two chunks with direct full-ensemble analysis.

    Assert that J, S, raw squeezing moments, and final squeezing are unchanged.
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
        n_processes=1,
    )
    expected = _compute_direct_moments(case)

    _assert_moments_allclose(actual, expected)


def test_combined_ensemble_preserves_trajectory_order(
    squeezing_post_analysis_case,
    squeezing_post_analysis_module,
):
    """
    Load uneven trajectory chunks in array-index order.

    Assert that every trajectory remains ordered and its seed records its chunk.
    """
    case = squeezing_post_analysis_case

    combined, _ = squeezing_post_analysis_module._load_combined_ensemble(
        case.input_dir,
        case.filename,
        len(case.chunk_sizes),
    )

    assert [
        trajectory.total_step_count for trajectory in combined.trajectories
    ] == [100, 101, 102]
    assert combined.seeds == [(0, 0), (1, 1), (1, 2)]
