# src/treetag/scoring.py
from __future__ import annotations
import pandas as pd
from .tree import init_tree

def subscores(
    root_cell: str,
    adata,
    tree_yaml: str = "tree.yaml",
    markers_yaml: str | None = None,  # None = faster; we only need structure
    only_leaves: bool = False,
    return_df: bool = False,
    require_all: bool = False,
):
    """
    Return the existing '<node>_score' columns for the subtree under `root_cell`.
    No computation—just collects columns already present in adata.obs.
    """
    # 1) Build subtree under root_cell
    G = init_tree(tree_yaml, markers_yaml=markers_yaml, root=root_cell)
    r = G.vs.find(name=root_cell).index
    desc = [i for i in G.subcomponent(r, mode="OUT") if i != r]
    if only_leaves:
        desc = [i for i in desc if G.outdegree(i) == 0]

    nodes = [G.vs[i]["name"] for i in desc]
    cols  = [f"{n}_score" for n in nodes]

    # 2) Keep only columns that actually exist
    existing = [c for c in cols if c in adata.obs]
    missing  = [c for c in cols if c not in adata.obs]

    if require_all and missing:
        raise KeyError(f"Missing expected score columns: {missing}")

    if return_df:
        return adata.obs[existing].copy() if existing else pd.DataFrame(index=adata.obs.index)

    return existing
