import os, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, treetag as tt

def tiny():
    # 9×6 toy: 3 families × 3 cells, 2 markers/family
    var = ["CD3D","TRAC","MS4A1","CD79A","LYZ","S100A8"]
    obs = [f"c{i}" for i in range(9)]
    X = np.array([
        [5,4,0,0,1,0],
        [6,3,0,0,0,1],
        [4,5,0,0,1,0],
        [0,0,6,5,0,1],
        [0,1,5,6,0,0],
        [1,0,4,5,0,1],
        [1,0,0,0,6,5],
        [0,1,0,0,5,6],
        [1,0,0,1,5,5],
    ], dtype=float)
    return ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=obs),
        var=pd.DataFrame(index=var),
    )

def _write_yaml(tmp_path):
    (tmp_path/"PBMC_tree.yaml").write_text(
        "root:\n  T:\n    T_Naive:\n  B:\n    B_Naive:\n  Myeloid:\n"
    )
    (tmp_path/"PBMC_markers.yaml").write_text(
        "T: [CD3D, TRAC]\nT_Naive: [CD3D]\n"
        "B: [MS4A1, CD79A]\nB_Naive: [MS4A1]\n"
        "Myeloid: [LYZ, S100A8]\n"
    )
    return str(tmp_path/"PBMC_tree.yaml"), str(tmp_path/"PBMC_markers.yaml")

def test_find_doublets_basic_runs(tmp_path):
    tree, markers = _write_yaml(tmp_path)
    a = tiny()

    # Provide the family-level scores expected by find_doublets
    X = a.X
    def col(g):
        j = a.var_names.get_loc(g)
        return X[:, j].toarray().ravel() if sp.issparse(X) else np.asarray(X[:, j]).ravel()
    a.obs["T_score"] = col("CD3D")
    a.obs["B_score"] = col("MS4A1")
    a.obs["Myeloid_score"] = col("LYZ")

    info = tt.find_doublets(a, tree_yaml=tree, markers_yaml=markers, root="root")

    assert {"doublet_score","cell#1","cell#2"} <= set(a.obs.columns)
    assert isinstance(info, dict)
