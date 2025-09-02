import igraph as ig
from ._init_tree import _init_tree as init_tree

def markers(cell_type, sign="pos", markers_yaml=None, adata=None):
    """
    markers_yaml format: cell_type: [gene1, gene2, -gene3]
    Build a single-node graph (cell_type), run init_tree(graph, markers_yaml),
    and return markers by sign: 'pos' | 'neg' | 'both'.
    """
    if sign not in {"pos", "neg", "both"}:
        raise ValueError("sign must be 'pos', 'neg', or 'both'")
    if not markers_yaml:
        raise ValueError("markers_yaml path required")

    # single-node ontology
    G = ig.Graph(directed=True)
    G.add_vertex(name=cell_type)

    # normalize via init_tree (splits into pos_markers / neg_markers, filters vs adata)
    G = init_tree(G, markers_yaml, root=cell_type, adata=adata)

    v = G.vs.find(name=cell_type)
    pos = list(v["pos_markers"] or [])
    neg = list(v["neg_markers"] or [])

    if adata is not None:
        var_names = adata.raw.var_names if getattr(adata, "raw", None) is not None else adata.var_names
        pos = [g for g in pos if g in var_names]
        neg = [g for g in neg if g in var_names]

    if sign == "both":
        return {"pos": pos, "neg": neg}
    return pos if sign == "pos" else neg
