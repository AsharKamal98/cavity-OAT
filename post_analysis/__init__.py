"""Post-simulation analysis package."""

from post_analysis.dressed_basis import build_dressed_basis_transformation
from post_analysis.generator_moments import compute_mean_generator_moments
from post_analysis.generator_rotation import (
    compute_gell_mann_rotation,
    rotate_generator_moments,
)
from post_analysis.harmonic_analysis import (
    compute_harmonic_analysis,
    compute_windowed_harmonic_analysis,
    extract_family_phase_series,
)
from post_analysis.j_moments import (
    J_FROM_GENERATORS,
    compute_j_components,
    compute_j_moments,
)
from post_analysis.j_modes import compute_j_modes, compute_j_modes_from_j_moments
from post_analysis.mfe_residuals import compute_mfe_residuals
from post_analysis.s_moments import compute_s_moments
from post_analysis.squeezing import compute_squeezing
from post_analysis.theory_benchmarks import phase1_ss_angles_for_nj
from post_analysis.uncertainty import (
    compute_sample_covariance,
    propagate_covariance,
    std_from_covariance,
)

__all__ = [
    "J_FROM_GENERATORS",
    "build_dressed_basis_transformation",
    "compute_gell_mann_rotation",
    "compute_harmonic_analysis",
    "compute_windowed_harmonic_analysis",
    "compute_mean_generator_moments",
    "compute_sample_covariance",
    "compute_j_moments",
    "compute_j_components",
    "compute_j_modes",
    "compute_j_modes_from_j_moments",
    "compute_mfe_residuals",
    "compute_s_moments",
    "compute_squeezing",
    "extract_family_phase_series",
    "phase1_ss_angles_for_nj",
    "propagate_covariance",
    "rotate_generator_moments",
    "std_from_covariance",
]
