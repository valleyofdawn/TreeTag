from typing import Literal, List, Optional
from ._init_tree import _init_tree as init_tree

def markers(cell_type: str,
            sign: str = "pos",
            markers_yaml: str = "markers.yaml",
            tree_yaml: str = "tree.yaml",
            adata=None):
    """
    Return the marker list for a given cell_type and sign ('pos' or 'neg'),
    using the new graph structure where markers are stored on nodes.
    If `adata` is provided, markers are filtered to genes present in the data.
    """
    if sign not in {"pos", "neg"}:
        raise ValueError("sign must be 'pos' or 'neg'")

    # Build graph with markers attached (filtered if adata is given)
    G = init_tree(tree_yaml, markers_yaml, root="root", adata=adata)  # use your actual root name

    # Look up the node
    try:
        v = G.vs.find(name=cell_type)
    except ValueError:
        print(f"Cell type '{cell_type}' not found in the tree. Skipping.")
        return []

    attr = "pos_markers" if sign == "pos" else "neg_markers"
    lst = list(v[attr]) if v[attr] is not None else []

    # Extra safety: if adata was provided, ensure presence in the matrix
    if adata is not None:
        var_names = adata.raw.var_names if getattr(adata, "raw", None) is not None else adata.var_names
        present = [g for g in lst if g in var_names]
        missing = [g for g in lst if g not in var_names]
        for g in missing:
            print(f"Marker '{g}' not found in data. Skipping.")
        return present

    return lst
