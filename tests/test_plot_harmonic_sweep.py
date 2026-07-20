import unittest

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from common.plotting.harmonic_analysis import plot_harmonic_sweep


class HarmonicSweepPlotTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_plots_frequency_and_thd_panels(self):
        parameter_values = np.array((0.0, 0.5, 1.0))
        frequencies = np.array((1.0, 1.2, 1.5))
        distortions = np.array((0.0, 0.1, 0.25))
        rms_amplitudes = np.array((0.2, 0.3, 0.4))
        offsets = np.array((0.5, 0.4, 0.3))

        fig, axes = plot_harmonic_sweep(
            parameter_values,
            frequencies,
            distortions,
            rms_amplitudes,
            offsets,
            labels=(r"$J_\perp$",),
            parameter_label=r"$\delta$",
        )

        self.assertEqual(axes.size, 4)
        np.testing.assert_allclose(axes[0].lines[0].get_xdata(), parameter_values)
        np.testing.assert_allclose(axes[0].lines[0].get_ydata(), frequencies)
        np.testing.assert_allclose(axes[1].lines[0].get_ydata(), rms_amplitudes)
        np.testing.assert_allclose(axes[2].lines[0].get_ydata(), offsets)
        np.testing.assert_allclose(axes[3].lines[0].get_ydata(), distortions)
        self.assertEqual(axes[0].get_ylabel(), r"$f_0$")
        self.assertEqual(axes[1].get_ylabel(), r"$A_{\mathrm{RMS}}$")
        self.assertEqual(axes[2].get_ylabel(), r"$c$")
        self.assertEqual(axes[3].get_ylabel(), "THD")
        self.assertEqual(axes[2].get_xlabel(), r"$\delta$")
        self.assertTrue(
            axes[0].get_shared_x_axes().joined(axes[0], axes[3])
        )
        self.assertEqual(fig.legends[0].get_texts()[0].get_text(), r"$J_\perp$")

    def test_overlays_on_supplied_axes(self):
        parameter_values = np.array((0.0, 1.0))
        fig, axes = plot_harmonic_sweep(
            parameter_values,
            (1.0, 1.1),
            (0.1, 0.2),
            (0.3, 0.4),
            (0.5, 0.6),
            labels=("bright",),
            parameter_label=r"$\delta$",
        )

        overlay_fig, overlay_axes = plot_harmonic_sweep(
            parameter_values,
            (0.8, 0.9),
            (0.2, 0.3),
            (0.4, 0.5),
            (0.6, 0.7),
            labels=("dark",),
            parameter_label=r"$\delta$",
            axes=axes,
            linestyle="--",
        )

        self.assertIs(overlay_fig, fig)
        self.assertEqual(len(overlay_axes[0].lines), 2)
        self.assertEqual(len(overlay_axes[1].lines), 2)
        self.assertEqual(len(overlay_axes[2].lines), 2)
        self.assertEqual(len(overlay_axes[3].lines), 2)
        self.assertEqual(
            overlay_axes[0].lines[0].get_color(),
            overlay_axes[0].lines[1].get_color(),
        )
        self.assertEqual(
            tuple(text.get_text() for text in fig.legends[0].get_texts()),
            ("bright", "dark"),
        )

    def test_rejects_invalid_sweep_arrays(self):
        invalid_inputs = (
            ((), (), (), "must not be empty"),
            ((0.0, 1.0), (1.0,), (0.1, 0.2), "matching lengths"),
            ((0.0, np.nan), (1.0, 1.1), (0.1, 0.2), "must be finite"),
            (((0.0, 1.0),), ((1.0, 1.1),), ((0.1, 0.2),), "one-dimensional"),
        )

        for parameters, frequencies, distortions, message in invalid_inputs:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    plot_harmonic_sweep(
                        parameters,
                        frequencies,
                        distortions,
                        distortions,
                        distortions,
                        labels=("curve",),
                        parameter_label="parameter",
                    )


if __name__ == "__main__":
    unittest.main()
