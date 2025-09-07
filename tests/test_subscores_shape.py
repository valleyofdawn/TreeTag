# tests/test_subscores_shape.py
import pytest, pandas as pd, numpy as np, treetag as tt


@pytest.mark.skipif(not hasattr(tt, "subscores"), reason="subscores API not present")
def test_subscores_output_type_and_shape(tiny_adata_csr, toy_yaml_pair):
    tree_y, mark_y = toy_yaml_pair

    out = tt.subscores(
        adata=tiny_adata_csr, tree_yaml=tree_y, markers_yaml=mark_y, root="root"
    )

    # Case 1: DataFrame
    if isinstance(out, pd.DataFrame):
        assert out.shape[0] in (0, tiny_adata_csr.n_obs)
        return

    # Case 2: list of dicts
    if isinstance(out, list) and out and isinstance(out[0], dict):
        df = pd.DataFrame(out)
        assert df.shape[0] in (0, tiny_adata_csr.n_obs)
        return

    # Case 3: list/array of numbers
    arr = np.asarray(out)
    assert arr.ndim >= 1
    # Accept either empty or length = n_obs
    assert arr.shape[0] in (0, tiny_adata_csr.n_obs)
