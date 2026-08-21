import numpy as np
import pytest

from solvers.mcwf.sim import _simulate_single_trajectory


def _assert_complete_moment_output(result, case):
    num_groups = len(case.Ni)
    generator_count = 8 * num_groups
    packed_count = generator_count * (generator_count + 1) // 2

    assert result.moments.num_trajectories == 1
    np.testing.assert_allclose(result.moments.t, case.t_eval)
    assert result.moments.first_moments.shape == (
        1,
        case.t_eval.size,
        num_groups,
        8,
    )
    assert result.moments.second_moments.shape == (
        1,
        case.t_eval.size,
        packed_count,
    )
    assert np.all(np.isfinite(result.moments.first_moments))
    assert np.all(np.isfinite(result.moments.second_moments))


@pytest.mark.benchmark(group="mcwf-single-trajectory")
def test_no_jump_trajectory(
    benchmark,
    no_jump_trajectory_case,
    no_jump_trajectory_precomputed,
):
    """
    Benchmark a fixed-seed trajectory dominated by precomputed full steps.

    Assert that no jump or partial-step fallback occurs and every requested
    moment row is returned.
    """
    case = no_jump_trajectory_case

    result = benchmark(
        _simulate_single_trajectory,
        Ni=case.Ni,
        omega_i=case.omega_i,
        Gamma=case.Gamma,
        integration_phases=case.integration_phases,
        sector_coeffs=case.sector_coeffs,
        dt=case.dt,
        t_eval=case.t_eval,
        seed_sequence=np.random.SeedSequence(case.seed),
        shifted_jump_operator=case.shifted_jump_operator,
        precomputed=no_jump_trajectory_precomputed,
    )

    _assert_complete_moment_output(result, case)
    assert result.diagnostics.jump_count == 0
    assert result.diagnostics.non_precomputed_step_count == 0

    benchmark.extra_info.update(
        atom_groups=case.Ni,
        num_sectors=len(case.sector_coeffs),
        num_saved_times=len(case.t_eval),
        jump_count=result.diagnostics.jump_count,
        total_step_count=result.diagnostics.total_step_count,
        non_precomputed_step_count=(
            result.diagnostics.non_precomputed_step_count
        ),
    )


@pytest.mark.benchmark(group="mcwf-single-trajectory")
def test_jump_heavy_trajectory(
    benchmark,
    jump_heavy_trajectory_case,
    jump_heavy_trajectory_precomputed,
):
    """
    Benchmark a fixed-seed trajectory with repeated jump localization.

    Assert that repeated jumps and partial-step propagations occur in the
    explicitly unshifted picture and every requested moment row is returned.
    """
    case = jump_heavy_trajectory_case

    result = benchmark(
        _simulate_single_trajectory,
        Ni=case.Ni,
        omega_i=case.omega_i,
        Gamma=case.Gamma,
        integration_phases=case.integration_phases,
        sector_coeffs=case.sector_coeffs,
        dt=case.dt,
        t_eval=case.t_eval,
        seed_sequence=np.random.SeedSequence(case.seed),
        shifted_jump_operator=case.shifted_jump_operator,
        precomputed=jump_heavy_trajectory_precomputed,
    )

    _assert_complete_moment_output(result, case)
    assert result.diagnostics.jump_count >= 20
    assert result.diagnostics.non_precomputed_step_count > 0

    benchmark.extra_info.update(
        atom_groups=case.Ni,
        num_sectors=len(case.sector_coeffs),
        num_saved_times=len(case.t_eval),
        jump_count=result.diagnostics.jump_count,
        total_step_count=result.diagnostics.total_step_count,
        non_precomputed_step_count=(
            result.diagnostics.non_precomputed_step_count
        ),
    )
