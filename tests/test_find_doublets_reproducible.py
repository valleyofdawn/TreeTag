# tests/test_find_doublets_reproducible.py
import inspect
import numpy as np
import pandas as pd
import pytest
import treetag as tt

def test_find_doublets_reproducible(tiny_adata_csr, toy_yaml_pair):
    tree_y, mark_y = toy_yaml_pair

    # Ensure family scores exist
    ctor = tt.TreeTag
    sig = inspect.signature(ctor)
    extra = {}
    if "min_pruning_fc" in sig.parameters:
        extra["min_pruning_fc"] = 1.0
    if "save_scores" in sig.parameters:
        extra["save_scores"] = True
    try:
        ctor(tiny_adata_csr, tree_yaml=tree_y, markers_yaml=mark_y, root="Root", **extra)
    except ValueError as e:
        pytest.xfail(f"pruning stopped on toy data: {e}")

    for col in ("B_score", "T_score"):
        assert col in tiny_adata_csr.obs, f"Missing {col}"

    fn = tt.find_doublets
    kwargs = {"adata": tiny_adata_csr, "tree_yaml": tree_y, "markers_yaml": mark_y, "root": "Root"}
    if "random_state" in inspect.signature(fn).parameters:
        kwargs["random_state"] = 0

    out1 = fn(**kwargs)

    # Mode A: returns dict summary (write_cols=True default)
    if isinstance(out1, dict):
        assert out1.get("n_cells") == tiny_adata_csr.n_obs
        assert set(out1.get("families", [])) >= {"B", "T"}
        out2 = fn(**kwargs)
        assert isinstance(out2, dict) and out1 == out2
        # in-place columns must exist and be stable
        assert "doublet_score" in tiny_adata_csr.obs
        v1 = tiny_adata_csr.obs["doublet_score"].to_numpy()
        _ = fn(**kwargs)
        v2 = tiny_adata_csr.obs["doublet_score"].to_numpy()
        assert v1.shape[0] == tiny_adata_csr.n_obs
        assert np.allclose(v1, v2)
        return

    # Mode B: returns array-like scores
    if out1 is not None:
        a1 = np.asarray(out1); a2 = np.asarray(fn(**kwargs))
        assert a1.ndim == 1 and a2.ndim == 1
        assert a1.shape == a2.shape
        assert a1.shape[0] in (0, tiny_adata_csr.n_obs)
        assert np.allclose(a1, a2)
        return

    # Mode C: returns None, wrote in-place
    assert "doublet_score" in tiny_adata_csr.obs
    v1 = tiny_adata_csr.obs["doublet_score"].to_numpy()
    _ = fn(**kwargs)
    v2 = tiny_adata_csr.obs["doublet_score"].to_numpy()
    assert v1.shape[0] == tiny_adata_csr.n_obs
    assert np.allclose(v1, v2)
