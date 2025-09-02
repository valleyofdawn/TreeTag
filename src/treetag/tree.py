def init_tree(
    tree_yaml: str | None,                  # now optional
    markers_yaml: str | None = None,        # already optional is fine
    root: str = "root",
    adata=None,
):
    """Subtree-aware loader (uses packaged PBMC YAMLs if paths are None)."""
    import yaml, igraph as ig
    from importlib.resources import files

    def _load_yaml(src, default_name=None):
        if src is None:
            if default_name is None:
                return {}
            res = files("treetag").joinpath("data", default_name)
            with res.open("r", encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        if hasattr(src, "open"):  # importlib.resources Traversable
            with src.open("r", encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        with open(src, "r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}

    # Use your PBMC files by default
    tree = _load_yaml(tree_yaml,     default_name="PBMC_tree.yaml")
    marker_dict = _load_yaml(markers_yaml, default_name="PBMC_markers.yaml")

    if not isinstance(tree, dict):
        raise ValueError("Tree YAML must be a nested mapping (dict of dicts).")

    # Flatten YAML -> edge list, skipping '!' branches
    def iter_edges(parent, mapping):
        for child, sub in (mapping or {}).items():
            if not child or str(child).startswith("_"):
                continue
            c = str(child)
            yield (parent, c)
            if isinstance(sub, dict):
                yield from iter_edges(c, sub)

    edges = []
    for top, sub in tree.items():
        if not top or str(top).startswith("_"):
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