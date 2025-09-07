# tests/test_treetag_bad_inputs.py
import pytest, treetag as tt


def test_treetag_missing_yaml_raises(tiny_adata_csr, toy_yaml_pair):
    _, mark_y = toy_yaml_pair
    with pytest.raises((FileNotFoundError, ValueError)):
        tt.TreeTag(
            tiny_adata_csr, tree_yaml="does_not_exist.yaml", markers_yaml=mark_y
        ).fit_predict()


def test_treetag_unknown_root_raises(tiny_adata_csr, toy_yaml_pair):
    tree_y, mark_y = toy_yaml_pair
    with pytest.raises((KeyError, ValueError)):
        tt.TreeTag(
            tiny_adata_csr, tree_yaml=tree_y, markers_yaml=mark_y, root="Nope"
        ).fit_predict()
