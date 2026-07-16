from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any

from parser.common import PhaseProtocol
from parser.j_moments import JMomentSeries


def save_j_moments_artifact(
    j_moments: JMomentSeries,
    phase_protocol: PhaseProtocol,
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump({"J": j_moments, "phase_protocol": phase_protocol}, fh)


def load_j_moments_artifact(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        return pickle.load(fh)
