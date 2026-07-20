import numpy as np
import pytest

from solvers.mcwf.sim import (
    propagate_blocks,
    propagate_blocks_with_propagators,
)


def _assert_matching_blocks(result, expected):
    assert len(result) == len(expected)
    assert all(output.shape == source.shape for output, source in zip(result, expected))
    assert all(np.all(np.isfinite(output)) for output in result)


@pytest.mark.benchmark(group="mcwf-propagation")
def test_full_step_propagation(benchmark, propagation_case):
    """
    Benchmark the hot full-step path using precomputed propagators.

    Assert that every sector produces a finite block with the expected shape.
    """
    case = propagation_case

    result = benchmark(
        propagate_blocks_with_propagators,
        case.psi_blocks,
        case.propagators,
    )

    _assert_matching_blocks(result, case.psi_blocks)
    benchmark.extra_info.update(
        num_sectors=len(case.psi_blocks),
        max_sector_dimension=max(block.size for block in case.psi_blocks),
    )


@pytest.mark.benchmark(group="mcwf-propagation")
def test_partial_step_propagation(benchmark, propagation_case):
    """
    Benchmark variable-step propagation through SciPy's exponential action.

    Assert that every sector produces a finite block with the expected shape.
    """
    case = propagation_case

    result = benchmark(
        propagate_blocks,
        case.psi_blocks,
        case.generators,
        case.partial_step,
    )

    _assert_matching_blocks(result, case.psi_blocks)
    benchmark.extra_info.update(
        num_sectors=len(case.psi_blocks),
        max_sector_dimension=max(block.size for block in case.psi_blocks),
        partial_step=case.partial_step,
    )
