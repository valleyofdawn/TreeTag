# src/treetag/doublets.py
import numpy as np
import pandas as pd
from ._init_tree import _init_tree as init_tree

def find_doublets(
    adata,
    tree_yaml: str,
    markers_yaml: str | None,
    root: str,
    write_cols: bool = True,
    eps: float = 0.1,
):
    """
    Compute simple doublet diagnostics from the first split (children of `root`).

    How it works
    ------------
    For each cell, look at the per-family scores written by TreeTag (one column per
    child: "<child>_score"). Let M1 be the highest score (winning family) and M2
    the second-highest. Define:
        doublet_score = (eps + M2) / (eps + M1)
    Because M2 ≤ M1, this ratio is in (0, 1]; values near 1.0 suggest the cell is
    pulled strongly by two families (doublet-like), while values near 0 indicate
    a clear single-family winner.

    Parameters
    ----------
    adata : AnnData
        Object containing the "<child>_score" columns from TreeTag (save_scores=True).
    tree_yaml : str
        Path to the ontology YAML (used to discover the root's children).
    markers_yaml : str | None
        Path to markers YAML, or None to skip loading markers (faster).
    root : str
        Node name to treat as the split root.
    write_cols : bool, default True
        If True, writes results to `adata.obs`.
    eps : float, default 0.1
        Small stabilizer to avoid division by zero and smooth tiny values.

    Writes (if write_cols=True)
    ---------------------------
    adata.obs['doublet_score'] : float
        (eps + M2) / (eps + M1), in (0, 1]; higher ⇒ more doublet-like.
    adata.obs['doublet_partner']  : category
        Runner-up family.

    Returns
    -------
    dict
        Summary with {'n_cells', 'families', 'root'}.
    """
    # 1) tree and direct children of root
    G = init_tree(tree_yaml, markers_yaml=markers_yaml, root=root)
    u = G.vs.find(name=root).index
    child_idxs = G.successors(u)
    if len(child_idxs) < 2:
        raise ValueError(f"Root '{root}' has fewer than 2 children.")
    families = [G.vs[i]["name"] for i in child_idxs]

    # 2) Score matrix from those children
    cols = [f"{c}_score" for c in families]
    missing = [c for c in cols if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing expected score columns from TreeTag: {missing}. Make sure TreeTag has been run with save_scores=True" )
    S = adata.obs[cols].to_numpy(dtype=float)
    S = np.nan_to_num(S, nan=0.0)

    # 3) Top-2 per cell
    n_cells = S.shape[0]
    idx_top2 = np.argpartition(S, -2, axis=1)[:, -2:]          # unsorted top-2 indices per row
    row = np.arange(n_cells)[:, None]
    top2_vals = S[row, idx_top2]                               # (n, 2)
    order = np.argsort(top2_vals, axis=1)                      # ascending within the pair
    best_idx   = idx_top2[row, order[:, 1][:, None]].ravel()   # index of M1
    second_idx = idx_top2[row, order[:, 0][:, None]].ravel()   # index of M2

    M1 = S[np.arange(n_cells), best_idx]
    M2 = S[np.arange(n_cells), second_idx]
    ratio = (eps + M2) / (eps + M1)
    top1 = np.array(families, dtype=object)[best_idx]
    top2 = np.array(families, dtype=object)[second_idx]

    if write_cols:
        adata.obs["doublet_score"] = ratio
        adata.obs["cell#1"]  = pd.Categorical(top1, categories=families)
        adata.obs["cell#2"]  = pd.Categorical(top2, categories=families)

    return {"n_cells": int(n_cells), "families": families, "root": root}
