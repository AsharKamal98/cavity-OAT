from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any, Sequence

from parser.common import Phase
from parser.j_moments import JMomentSeries


def save_j_moments_artifact(
    j_moments: JMomentSeries,
    phases: Sequence[Phase],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump({"J": j_moments, "phases": list(phases)}, fh)


def load_j_moments_artifact(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        return pickle.load(fh)
