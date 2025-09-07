import os
import tempfile
import yaml
from ._init_tree import _init_tree as init_tree


def markers(cell_type, sign="pos", markers_yaml=None, adata=None):
    """
    markers_yaml format: cell_type: [gene1, gene2, -gene3]
    Build a one-node ontology YAML, run init_tree, and return markers.
    """
    if sign not in {"pos", "neg", "both"}:
        raise ValueError("sign must be 'pos', 'neg', or 'both'")
    if not markers_yaml:
        raise ValueError("markers_yaml path required")

    # build one-node ontology yaml
    tiny_tree = {cell_type: None}
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as ft:
        yaml.safe_dump(tiny_tree, ft)
        tree_path = ft.name

    try:
        G = init_tree(tree_path, markers_yaml, root=cell_type, adata=adata)
        v = G.vs.find(name=cell_type)
        pos = list(v["pos_markers"] or [])
        neg = list(v["neg_markers"] or [])
        if adata is not None:
            var_names = (
                adata.raw.var_names
                if getattr(adata, "raw", None) is not None
                else adata.var_names
            )
            pos = [g for g in pos if g in var_names]
            neg = [g for g in neg if g in var_names]
        if sign == "both":
            return {"pos": pos, "neg": neg}
        return pos if sign == "pos" else neg
    finally:
        try:
            os.unlink(tree_path)
        except Exception:
            pass
