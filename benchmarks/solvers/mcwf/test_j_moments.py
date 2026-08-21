import pytest

from post_analysis import compute_j_moments, compute_mean_generator_moments


def _compute_mean_and_j_moments(case):
    mean_moments = compute_mean_generator_moments(case.moments)
    return compute_j_moments(mean_moments, metadata=case.metadata)


@pytest.mark.benchmark(group="mcwf-moment-post-analysis")
def test_mean_and_j_moments_many_snapshots(benchmark, moment_analysis_case):
    """
    Benchmark trajectory averaging and generic J conversion on a large grid.

    Assert that the returned aggregated and group moment series contain one
    value for every saved snapshot.
    """
    case = moment_analysis_case

    moments = benchmark(_compute_mean_and_j_moments, case)

    num_snapshots = case.moments.t.size
    assert len(moments.t) == num_snapshots
    assert moments.x is not None
    assert all(len(group) == num_snapshots for group in moments.x)

    benchmark.extra_info.update(
        num_trajectories=case.moments.num_trajectories,
        num_snapshots=num_snapshots,
        num_groups=case.moments.first_moments.shape[2],
    )
