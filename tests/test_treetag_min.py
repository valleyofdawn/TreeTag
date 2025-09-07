# tests/test_treetag_min.py
import pandas as pd, scipy.sparse as sp, anndata as ad, treetag as tt


def test_public_api():
    assert hasattr(tt, "TreeTag") and hasattr(tt, "markers")


def test_markers_returns_genes(yaml_paths):
    _, markers = yaml_paths
    out = tt.markers("B", markers_yaml=markers)
    assert isinstance(out, list) and "MS4A1" in out


def test_treetag_labels_inplace_or_return(tiny_ad, yaml_paths):
    tree, markers = yaml_paths
    a = tiny_ad.copy()
    res = tt.TreeTag(
        adata=a, tree_yaml=tree, markers_yaml=markers, root="root", verbose=False
    )
    b = res if hasattr(res, "obs") else a
    assert "TreeTag" in b.obs and b.obs["TreeTag"].notna().all()
