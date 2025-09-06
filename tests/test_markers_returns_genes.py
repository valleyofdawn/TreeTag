import pandas as pd, numpy as np, anndata as ad, scipy.sparse as sp, treetag as tt

def test_markers_returns_genes(yaml_paths):
    _, markers = yaml_paths
    out = tt.markers("B", markers_yaml=markers)
    assert "MS4A1" in out