import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, pytest

@pytest.fixture
def tiny_ad():
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
    ], float)
    return ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=obs), var=pd.DataFrame(index=var))

@pytest.fixture
def yaml_paths(tmp_path):
    tree = tmp_path/"PBMC_tree.yaml"
    markers = tmp_path/"PBMC_markers.yaml"
    tree.write_text("root:\n  T:\n    T_Naive:\n  B:\n    B_Naive:\n  Myeloid:\n")
    markers.write_text(
        "T: [CD3D, TRAC]\nT_Naive: [CD3D]\n"
        "B: [MS4A1, CD79A]\nB_Naive: [MS4A1]\n"
        "Myeloid: [LYZ, S100A8]\n"
    )
    return str(tree), str(markers)
