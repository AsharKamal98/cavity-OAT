import pytest

from solvers.mcwf.sim import build_precomputed_trajectory_data


@pytest.mark.slow
@pytest.mark.benchmark(group="mcwf-precompute")
def test_precompute_two_group(benchmark, precompute_case):
    """
    Benchmark construction of reusable two-group MCWF data.

    Assert that the result contains one aligned propagator for every
    integration phase and populated sector.
    """
    case = precompute_case

    precomputed = benchmark(
        build_precomputed_trajectory_data,
        Ni=case.Ni,
        omega_i=case.omega_i,
        Gamma=case.Gamma,
        integration_phases=case.integration_phases,
        sector_coeffs=case.sector_coeffs,
        dt=case.dt,
    )

    num_sectors = len(case.sector_coeffs)
    assert len(precomputed["sector_list"]) == num_sectors
    assert len(precomputed["integration_phase_propagators"]) == len(
        case.integration_phases
    )
    assert all(
        len(propagators) == num_sectors
        for propagators in precomputed["integration_phase_propagators"]
    )

    benchmark.extra_info.update(
        atom_groups=case.Ni,
        num_sectors=num_sectors,
        max_sector_dimension=max(precomputed["dims"].values()),
        num_integration_phases=len(case.integration_phases),
    )
