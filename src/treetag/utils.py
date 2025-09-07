# src/treetag/utils.py
from importlib.resources import files
from pathlib import Path
from collections.abc import Sequence


def list_files() -> list[str]:
    d = files("treetag") / "data"
    return sorted(p.name for p in d.iterdir() if p.is_file())


def fetch_files(
    names: str | Path | Sequence[str | Path] | None = None,
    dest: str | Path = ".",
    overwrite: bool = False,
) -> list[str]:
    """Copy 1 or many files from treetag/data to dest. Returns list of dest paths."""
    if names is None:
        names = list_files()
    elif isinstance(names, (str, Path)):
        names = [names]
    # normalize to bare filenames
    names = [Path(n).name for n in names]

    available = set(list_files())
    missing = [n for n in names if n not in available]
    if missing:
        raise ValueError(f"Not found: {missing}. Available: {sorted(available)}")

    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    out_paths = []
    for n in names:
        src = files("treetag") / "data" / n
        out = dest / n
        if overwrite or not out.exists():
            out.write_bytes(src.read_bytes())
        out_paths.append(str(out))
    return out_paths
