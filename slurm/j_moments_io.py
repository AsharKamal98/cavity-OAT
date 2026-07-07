from __future__ import annotations

import pickle
from pathlib import Path

from parser.j_moments import JMomentSeries


def save_j_moments_pickle(j_moments: JMomentSeries, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump(j_moments, fh)


def load_j_moments_pickle(path: Path) -> JMomentSeries:
    with path.open("rb") as fh:
        return pickle.load(fh)
