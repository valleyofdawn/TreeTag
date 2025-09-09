import pandas as pd, numpy as np, anndata as ad, scipy.sparse as sp, treetag as tt


def test_find_doublets_basic_runs(tiny_ad, yaml_paths):
    a = tiny_ad.copy()  # AnnData is mutable
    tree, markers = yaml_paths
    X = a.X

    def col(g):
        j = a.var_names.get_loc(g)
        return X[:, j].toarray().ravel()

    a.obs["T_score"] = col("CD3D")
    a.obs["B_score"] = col("MS4A1")
    a.obs["Myeloid_score"] = col("LYZ")
    info = tt.find_doublets(a, tree_yaml=tree, markers_yaml=markers, root="root")
    assert {"mismatch_score", "cell#1", "cell#2"} <= set(a.obs.columns)
    assert isinstance(info, dict)
