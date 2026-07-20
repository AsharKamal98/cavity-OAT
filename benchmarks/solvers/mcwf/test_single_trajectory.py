import numpy as np
import pytest

from solvers.mcwf.sim import _simulate_single_trajectory, total_norm2_list


@pytest.mark.benchmark(group="mcwf-single-trajectory")
def test_no_jump_trajectory(
    benchmark,
    no_jump_trajectory_case,
    no_jump_trajectory_precomputed,
):
    """
    Benchmark a fixed-seed trajectory dominated by precomputed full steps.

    Assert that no jump or partial-step fallback occurs and the final state is
    normalized.
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

    np.testing.assert_allclose(
        [snapshot.time for snapshot in result.snapshots],
        case.t_eval,
    )
    assert result.jump_count == 0
    assert result.non_precomputed_step_count == 0
    assert np.isclose(total_norm2_list(result.final_sector_blocks.values()), 1.0)

    benchmark.extra_info.update(
        atom_groups=case.Ni,
        num_sectors=len(case.sector_coeffs),
        num_saved_times=len(case.t_eval),
        jump_count=result.jump_count,
        total_step_count=result.total_step_count,
        non_precomputed_step_count=result.non_precomputed_step_count,
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
    explicitly unshifted picture and the final state remains normalized.
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

    np.testing.assert_allclose(
        [snapshot.time for snapshot in result.snapshots],
        case.t_eval,
    )
    assert result.jump_count >= 20
    assert result.non_precomputed_step_count > 0
    assert np.isclose(total_norm2_list(result.final_sector_blocks.values()), 1.0)

    benchmark.extra_info.update(
        atom_groups=case.Ni,
        num_sectors=len(case.sector_coeffs),
        num_saved_times=len(case.t_eval),
        jump_count=result.jump_count,
        total_step_count=result.total_step_count,
        non_precomputed_step_count=result.non_precomputed_step_count,
    )
