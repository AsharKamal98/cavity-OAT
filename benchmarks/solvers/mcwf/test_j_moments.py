import pytest

from solvers.mcwf.j_moments import compute_mcwf_j_moments


@pytest.mark.slow
@pytest.mark.benchmark(group="mcwf-moment-extraction")
def test_moment_extraction_many_snapshots(benchmark, moment_extraction_ensemble):
    """
    Benchmark serial J-moment extraction on a large saved-time grid.

    Assert that the returned total and group-resolved moment series contain one
    value for every saved snapshot.
    """
    ensemble = moment_extraction_ensemble

    moments = benchmark(
        compute_mcwf_j_moments,
        ensemble,
        n_processes=1,
    )

    num_snapshots = len(ensemble.metadata.t_eval)
    assert len(moments.t) == num_snapshots
    assert moments.x_groups is not None
    assert all(len(group) == num_snapshots for group in moments.x_groups)

    benchmark.extra_info.update(
        num_trajectories=len(ensemble.trajectories),
        num_snapshots=num_snapshots,
        num_sectors=len(ensemble.metadata.sectors),
    )
