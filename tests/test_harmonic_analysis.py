import unittest

import numpy as np

from post_analysis.harmonic_analysis import compute_harmonic_analysis


def _phase_difference(actual, expected):
    return np.angle(np.exp(1j * (actual - expected)))


class HarmonicAnalysisTests(unittest.TestCase):
    def test_recovers_bin_aligned_harmonics(self):
        t = np.arange(0.0, 10.0, 0.01)
        values = (
            1.3
            + 2.0 * np.cos(2.0 * np.pi * 0.5 * t + 0.3)
            + 0.4 * np.cos(2.0 * np.pi * 1.0 * t - 0.7)
        )

        signal = compute_harmonic_analysis(
            t,
            values,
            max_harmonic=3,
        ).signals[0]

        self.assertAlmostEqual(signal.fundamental_frequency, 0.5, places=7)
        np.testing.assert_allclose(signal.amplitudes[:2], (2.0, 0.4), atol=1e-7)
        np.testing.assert_allclose(
            _phase_difference(signal.phase_offsets[:2], (0.3, -0.7)),
            0.0,
            atol=1e-7,
        )
        self.assertAlmostEqual(signal.total_harmonic_distortion, 0.2, places=7)

    def test_recovers_harmonics_with_a_partial_final_cycle(self):
        t = np.arange(0.0, 8.3, 0.01)
        values = (
            -0.2
            + 1.5 * np.cos(2.0 * np.pi * 0.7 * t - 0.4)
            + 0.3 * np.cos(2.0 * np.pi * 1.4 * t + 0.8)
        )

        signal = compute_harmonic_analysis(t, values, max_harmonic=2).signals[0]

        self.assertAlmostEqual(signal.fundamental_frequency, 0.7, places=7)
        np.testing.assert_allclose(signal.amplitudes, (1.5, 0.3), atol=1e-7)
        np.testing.assert_allclose(
            _phase_difference(signal.phase_offsets, (-0.4, 0.8)),
            0.0,
            atol=1e-7,
        )

    def test_accepts_multiple_series(self):
        t = np.arange(0.0, 6.0, 0.01)
        first = np.cos(2.0 * np.pi * 1.0 * t)
        second = 0.5 * np.cos(2.0 * np.pi * 1.5 * t + 0.2)

        result = compute_harmonic_analysis(
            t,
            (first, second),
            max_harmonic=2,
        )

        self.assertFalse(hasattr(result, "t"))
        np.testing.assert_allclose(
            tuple(signal.fundamental_frequency for signal in result.signals),
            (1.0, 1.5),
            atol=1e-7,
        )

    def test_rejects_nonuniform_time_grid(self):
        t = np.array((0.0, 0.1, 0.21, 0.3))
        with self.assertRaisesRegex(ValueError, "uniformly sampled"):
            compute_harmonic_analysis(t, np.ones_like(t))

    def test_returns_none_for_less_than_one_observed_oscillation(self):
        t = np.arange(0.0, 2.0, 0.01)
        values = np.cos(2.0 * np.pi * 0.2 * t)

        result = compute_harmonic_analysis(t, values)

        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
