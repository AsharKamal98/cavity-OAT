from __future__ import annotations

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from slurm.trajectory_ensemble_io import (
    load_combined_trajectory_ensemble,
    save_trajectory_ensemble_artifact,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine raw MCWF array chunks before post-analysis."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--filename", default="ensemble")
    parser.add_argument("--num-files", type=int, default=4)
    args = parser.parse_args()

    ensemble, metadata = load_combined_trajectory_ensemble(
        args.input_dir, args.filename, args.num_files,
    )
    output_path = args.input_dir / f"{args.filename}.pkl"
    temporary_path = output_path.with_suffix(".pkl.tmp")
    save_trajectory_ensemble_artifact(
        ensemble, metadata, temporary_path, array_index=-1,
    )
    temporary_path.replace(output_path)
    print(
        f"Combined {args.num_files} chunks into {output_path} "
        f"({ensemble.moments.num_trajectories} trajectories)."
    )


if __name__ == "__main__":
    main()
