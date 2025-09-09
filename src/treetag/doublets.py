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
    beta: float | None = None,
):
    """
    Compute simple doublet diagnostics from the first split (children of `root`).

    Method
    ------
    For each cell, read per-family scores ("<child>_score"). Let M1 be the max score
    and M2 the second max. Define M = M1 + M2 and:
        doublet_score = (M2 / (M1 + 1e-8)) * ( M / (M + beta) )
    with beta = median(M) by default. This preserves the M2/M1 ratio at normal signal
    and suppresses it when both M1 and M2 are tiny. No hard gating.

    Parameters
    ----------
    adata : AnnData
        Contains the "<child>_score" columns from TreeTag (save_scores=True).
    tree_yaml : str
        Path to ontology YAML (used to discover root's children).
    markers_yaml : str | None
        Path to markers YAML, or None to skip loading markers.
    root : str
        Node name to treat as the split root.
    write_cols : bool, default True
        If True, writes results to `adata.obs`.
    beta : float | None
        Shrinkage scale for low-signal suppression. If None, uses median(M).

    Writes (if write_cols=True)
    ---------------------------
    adata.obs['doublet_score'] : float
        (M2/(M1+1e-8)) * (M/(M+beta)); higher ⇒ more doublet-like.
    adata.obs['doublet_partner'] : category
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
        raise KeyError(
            f"Missing expected score columns from TreeTag: {missing}. "
            f"Run TreeTag with save_scores=True."
        )
    S = adata.obs[cols].to_numpy(dtype=float)
    S = np.nan_to_num(S, nan=0.0)

    # 3) Top-2 per cell
    n_cells = S.shape[0]
    idx_top2 = np.argpartition(S, -2, axis=1)[:, -2:]
    row = np.arange(n_cells)[:, None]
    top2_vals = S[row, idx_top2]  # (n, 2)
    order = np.argsort(top2_vals, axis=1)  # ascending within the pair
    best_idx = idx_top2[row, order[:, 1][:, None]].ravel()  # index of M1
    second_idx = idx_top2[row, order[:, 0][:, None]].ravel()  # index of M2

    M1 = S[np.arange(n_cells), best_idx]
    M2 = S[np.arange(n_cells), second_idx]
    M = M1 + M2

    # beta: scale for magnitude shrinkage
    if beta is None:
        med = np.median(M) if n_cells > 0 else 0.0
        beta = med if med > 0 else 1.0  # avoid zero beta

    ratio = M2 / (M1 + 1e-8)
    mag = M / (M + beta)
    score = ratio * mag

    top1 = np.array(families, dtype=object)[best_idx]
    top2 = np.array(families, dtype=object)[second_idx]

    if write_cols:
        adata.obs["mismatch_score"] = score
        adata.obs["cell#1"] = pd.Categorical(top1, categories=families)
        adata.obs["cell#2"] = pd.Categorical(top2, categories=families)

    return {"n_cells": int(n_cells), "families": families, "root": root}
