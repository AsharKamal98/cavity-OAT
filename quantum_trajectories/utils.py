import multiprocessing as mp
from typing import Iterable, Optional

from tqdm.auto import tqdm

def map_with_optional_pool(
    worker,
    items: Iterable,
    *,
    n_processes: Optional[int],
    progress_desc: str,
):
    """
    Run an independent worker over items serially or with multiprocessing.
    """
    items = list(items)
    if n_processes is None or n_processes == 1:
        return [worker(item) for item in tqdm(items, desc=progress_desc)]

    if n_processes == -1:
        n_processes = mp.cpu_count()
    if n_processes <= 0:
        raise ValueError("n_processes must be None, 1, -1, or a positive integer.")

    ctx = mp.get_context()
    with ctx.Pool(processes=n_processes) as pool:
        return list(tqdm(pool.imap(worker, items), total=len(items), desc=progress_desc))


__all__ = [
    "map_with_optional_pool",
]
