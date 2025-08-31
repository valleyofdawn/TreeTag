def init_tree(
    tree_yaml: str,                    # Path to ontology tree YAML (nested dict of nodes)
    markers_yaml: str | None = None,   # Optional markers YAML; if None, skip loading
    root: str = "root",                # Node name to treat as subtree root
    adata=None,                        # Optional AnnData; filters markers to dataset genes
):
    """Subtree-aware loader.

    - YAML may encode a full tree with any top-level key; `root` can be any node name inside it.
    - Builds the directed graph from the entire YAML, then slices to the subtree rooted at `root`.
    - If `markers_yaml` is provided, attaches v['pos_markers'], v['neg_markers'] (optionally filtered to `adata` genes).
      Marker YAML format: node -> [genes], with negatives prefixed by '-'.
    - If `markers_yaml` is None, no marker attributes are attached (faster for tasks like plot_tree).
    """
    import yaml, igraph as ig

    # Load tree
    with open(tree_yaml, "r") as f:
        tree = yaml.safe_load(f) or {}
    if not isinstance(tree, dict):
        raise ValueError("Tree YAML must be a nested mapping (dict of dicts).")

    # Flatten YAML -> edge list, skipping '!' branches
    def iter_edges(parent, mapping):
        for child, sub in (mapping or {}).items():
            if not child or str(child).startswith("!"):
                continue
            c = str(child)
            yield (parent, c)
            if isinstance(sub, dict):
                yield from iter_edges(c, sub)

    edges = []
    for top, sub in tree.items():
        if not top or str(top).startswith("!"):
            continue
        t = str(top)
        edges.extend(iter_edges(t, sub if isinstance(sub, dict) else {}))

    # Build full graph
    G_full = ig.Graph.TupleList(edges, directed=True, vertex_name_attr="name")
    # ensure isolated top-level vertices exist
    for top in tree.keys():
        nm = str(top)
        if not nm.startswith("_") and nm not in G_full.vs["name"]:
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
    with open(markers_yaml, "r") as f:
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
        p, q = to_pos_neg(marker_dict.get(n, []))
        if vg is not None:
            p = [g for g in p if g in vg]
            q = [g for g in q if g in vg]
        pos_attr.append(tuple(p))
        neg_attr.append(tuple(q))

    G.vs["pos_markers"] = pos_attr
    G.vs["neg_markers"] = neg_attr
    all_mark = [*pos_attr, *neg_attr]
    G["marker_union"] = tuple(sorted({g for lst in all_mark for g in lst}))
    return G
