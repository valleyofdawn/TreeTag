# tests/test_treetag_output_contract.py
import pandas as pd, treetag as tt


def test_treetag_returns_series_when_not_inplace(tiny_adata_csr, toy_yaml_pair):
    tree_y, mark_y = toy_yaml_pair

    # Call the public API directly; it may return None and write into obs
    res = tt.TreeTag(tiny_adata_csr, tree_yaml=tree_y, markers_yaml=mark_y, root="Root")

    if isinstance(res, pd.Series):
        ser = res
    else:
        # Expect in-place write
        assert "TreeTag" in tiny_adata_csr.obs
        ser = tiny_adata_csr.obs["TreeTag"]

    assert isinstance(ser, pd.Series)
    assert len(ser) == tiny_adata_csr.n_obs
    assert ser.index.equals(tiny_adata_csr.obs_names)
