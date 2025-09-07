import pytest, numpy as np, pandas as pd
import scipy.sparse as sp
import anndata as ad


@pytest.fixture(scope="function")
def tiny_adata_csr():
    # 9 cells × 6 genes: T, B, and myeloid blocks
    genes = ["CD3D", "TRAC", "MS4A1", "CD79A", "LYZ", "S100A8"]
    X = np.array(
        [
            [5, 4, 0, 0, 1, 0],  # T-like
            [6, 3, 0, 0, 0, 1],  # T-like
            [4, 5, 0, 0, 1, 0],  # T-like
            [0, 0, 6, 5, 0, 1],  # B-like
            [0, 1, 5, 6, 0, 0],  # B-like
            [1, 0, 4, 5, 0, 1],  # B-like
            [1, 0, 0, 0, 6, 5],  # Myeloid
            [0, 1, 0, 0, 5, 6],  # Myeloid
            [1, 0, 0, 1, 5, 5],  # Myeloid
        ],
        dtype=float,
    )
    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(X.shape[0])]),
        var=pd.DataFrame(index=genes),
    )
    adata.raw = adata
    return adata


@pytest.fixture(scope="function")
def toy_yaml_pair(tmp_path_factory):
    d = tmp_path_factory.mktemp("yaml")
    tree_p = d / "PBMC_tree.yaml"
    mark_p = d / "PBMC_markers.yaml"
    tree_p.write_text("Root:\n  B:\n  T:\n")
    # ≥2 +markers per family; explicit negatives vs. others
    mark_p.write_text(
        "B: [MS4A1, CD79A, -CD3D, -TRAC, -LYZ, -S100A8]\n"
        "T: [CD3D, TRAC, -MS4A1, -CD79A, -LYZ, -S100A8]\n"
    )
    return str(tree_p), str(mark_p)
