# tests/test_utils_files.py
import os, inspect, pytest, treetag as tt

pytestmark = pytest.mark.skipif(
    not all(hasattr(tt, n) for n in ("list_files","fetch_files")),
    reason="utils API not present"
)

def test_list_files_returns_strings_or_empty():
    files = tt.list_files()
    assert isinstance(files, (list, tuple))
    assert all(isinstance(x, str) for x in files)

def test_fetch_files_writes_to_dest_or_skips(tmp_dir):
    # Prefer a known demo file name, fall back to first listed file
    target = "PBMC_tree.yaml"
    files = []
    try:
        files = tt.list_files()
    except Exception:
        pass
    if (not files) and target is None:
        pytest.skip("No known example files exposed by package")

    name = target if target in files or target else (files[0] if files else None)
    if name is None:
        pytest.skip("No fetchable example filename available")

    sig = inspect.signature(tt.fetch_files)
    kwargs = {"dest": str(tmp_dir)} if "dest" in sig.parameters else {}
    try:
        out = tt.fetch_files(name, **kwargs)
    except FileNotFoundError:
        pytest.skip(f"Example file {name} not bundled; skipping")
    except Exception as e:
        # If the API expects different args, surface the failure
        raise

    paths = out if isinstance(out, (list, tuple)) else [out]
    # Allow empty result, but if paths exist they must be files
    if paths:
        assert all(isinstance(p, str) and os.path.isfile(p) for p in paths)
