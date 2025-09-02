def init_tree(
    tree_yaml: str,
    markers_yaml: str | None = None,
    root: str = "root",
    adata=None,
):
    """Subtree-aware loader.

    - YAML may encode a full tree with any top-level key; `root` can be any node name inside it.
    - Builds the directed graph from the entire YAML, then slices to the subtree at `root`.
    - If `markers_yaml` is provided, attaches v['pos_markers'], v['neg_markers'] (optionally gene-filtered).
      Marker YAML: node -> [genes], negatives prefixed by '-'. (Quote them in YAML.)
    - Nodes starting with '_' (or '!') are ignored/disabled.
    """
    import yaml, igraph as ig
    from pathlib import Path

    # Load tree
    p = Path(tree_yaml)
    with p.open("r", encoding="utf-8") as f:
        tree = yaml.safe_load(f) or {}
    if not isinstance(tree, dict):
        raise ValueError("Tree YAML must be a nested mapping (dict of dicts).")

    def _disabled(name: str) -> bool:
        s = str(name)
        return s.startswith("_")

    # Flatten YAML -> edge list, skipping disabled branches
    def iter_edges(parent, mapping):
        for child, sub in (mapping or {}).items():
            if not child or _disabled(child):
                continue
            c = str(child)
            yield (parent, c)
            if isinstance(sub, dict):
                yield from iter_edges(c, sub)

    edges = []
    for top, sub in tree.items():
        if not top or _disabled(top):
            continue
        t = str(top)
        edges.extend(iter_edges(t, sub if isinstance(sub, dict) else {}))

    # Build full graph
    G_full = ig.Graph.TupleList(edges, directed=True, vertex_name_attr="name")
    for top in tree.keys():
        nm = str(top)
        if (not _disabled(nm)) and nm not in G_full.vs["name"]:
            G_full.add_vertex(name=nm)

    if root not in G_full.vs["name"]:
        raise ValueError(f"Requested root '{root}' not found in YAML tree.")

    # Subtree via igraph
    r = G_full.vs.find(name=root).index
    idxs = G_full.subcomponent(r, mode="OUT")
    G = G_full.induced_subgraph(idxs)
    G["root"] = root

    # Early exit if no markers requested
    if markers_yaml is None:
        return G

    # Load markers (node -> list[str], negatives prefixed by '-')
    p = Path(markers_yaml)
    with p.open("r", encoding="utf-8") as f:
        marker_dict = yaml.safe_load(f) or {}
    if not isinstance(marker_dict, dict):
        marker_dict = {}

    # Optional gene filter
    vg = None
    if adata is not None:
        try:
            vg = set(adata.raw.var_names if getattr(adata, "raw", None) is not None else adata.var_names)
        except Exception:
            vg = None

    def to_pos_neg(lst):
        lst = list(lst or [])
        neg, pos = [], []
        for g in lst:
            s = str(g)
            if s.startswith("-"):
                neg.append(s[1:])
            else:
                pos.append(s)
        return pos, neg

    pos_attr, neg_attr = [], []
    for n in G.vs["name"]:
        p_genes, n_genes = to_pos_neg(marker_dict.get(n, []))
        if vg is not None:
            p_genes = [g for g in p_genes if g in vg]
            n_genes = [g for g in n_genes if g in vg]
        pos_attr.append(tuple(p_genes))
        neg_attr.append(tuple(n_genes))

    G.vs["pos_markers"] = pos_attr
    G.vs["neg_markers"] = neg_attr
    all_mark = [*pos_attr, *neg_attr]
    G["marker_union"] = tuple(sorted({g for lst in all_mark for g in lst}))
    return G