import unittest

import numpy as np

from parser.bloch_vector import BlochVectorSeries
from post_analysis.j_modes import compute_j_modes


class BlochVectorSeriesTests(unittest.TestCase):
    def test_derives_length_and_xy_length(self):
        vector = BlochVectorSeries(
            x=np.array((3.0, 0.0)),
            y=np.array((4.0, 0.0)),
            z=np.array((0.0, 5.0)),
        )

        np.testing.assert_allclose(vector.length, (5.0, 5.0))
        np.testing.assert_allclose(vector.xy_length, (5.0, 0.0))

    def test_j_modes_receive_derived_lengths(self):
        modes = compute_j_modes(
            np.array((0.0, 1.0)),
            np.array((3.0, 0.0)),
            np.array((4.0, 0.0)),
            np.array((0.0, 5.0)),
            populations=(1.0,),
            omega_groups=(1.0,),
        )

        np.testing.assert_allclose(modes.bright.length, (5.0, 5.0))
        np.testing.assert_allclose(modes.bright.xy_length, (5.0, 0.0))

    def test_rejects_mismatched_component_shapes(self):
        with self.assertRaisesRegex(ValueError, "matching one-dimensional"):
            BlochVectorSeries(
                x=np.array((1.0, 2.0)),
                y=np.array((1.0,)),
                z=np.array((1.0, 2.0)),
            )


if __name__ == "__main__":
    unittest.main()
