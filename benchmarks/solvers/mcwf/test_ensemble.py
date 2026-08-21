import numpy as np
import pytest

from solvers.mcwf.ensamble_sim import run_trajectory_ensemble


def _assert_complete_ensemble(result, case):
    num_groups = len(case.parameters.Ni)
    generator_count = 8 * num_groups
    packed_count = generator_count * (generator_count + 1) // 2

    assert result.moments.num_trajectories == case.ntraj
    assert len(result.seeds) == case.ntraj
    assert result.metadata is not None
    np.testing.assert_allclose(result.moments.t, case.t_eval)
    np.testing.assert_allclose(result.metadata.t_eval, case.t_eval)
    assert result.moments.first_moments.shape == (
        case.ntraj,
        case.t_eval.size,
        num_groups,
        8,
    )
    assert result.moments.second_moments.shape == (
        case.ntraj,
        case.t_eval.size,
        packed_count,
    )


def _add_ensemble_metadata(benchmark, result, case):
    benchmark.extra_info.update(
        atom_groups=case.parameters.Ni,
        num_trajectories=case.ntraj,
        num_processes=case.n_processes,
        chunksize=case.chunksize,
        num_saved_times=result.moments.t.size,
    )


@pytest.mark.benchmark(group="mcwf-ensemble")
def test_serial_ensemble(benchmark, serial_ensemble_case):
    """
    Benchmark the public end-to-end serial ensemble runner.

    Assert that every requested trajectory, seed, and saved-time grid is
    returned.
    """
    case = serial_ensemble_case

    result = benchmark(
        run_trajectory_ensemble,
        case.parameters,
        t_eval=case.t_eval,
        seed=case.seed,
        ntraj=case.ntraj,
        n_processes=case.n_processes,
        chunksize=case.chunksize,
    )

    _assert_complete_ensemble(result, case)
    _add_ensemble_metadata(benchmark, result, case)


@pytest.mark.expensive
@pytest.mark.benchmark(group="mcwf-ensemble")
def test_multiprocessing_ensemble(benchmark, multiprocessing_ensemble_case):
    """
    Benchmark end-to-end ensemble throughput with two worker processes.

    Assert that multiprocessing returns every requested trajectory, seed, and
    saved-time grid.
    """
    case = multiprocessing_ensemble_case

    result = benchmark.pedantic(
        run_trajectory_ensemble,
        args=(case.parameters,),
        kwargs={
            "t_eval": case.t_eval,
            "seed": case.seed,
            "ntraj": case.ntraj,
            "n_processes": case.n_processes,
            "chunksize": case.chunksize,
        },
        rounds=3,
        iterations=1,
    )

    _assert_complete_ensemble(result, case)
    _add_ensemble_metadata(benchmark, result, case)
