from .tree import init_tree
from .convert import convert

def TreeTag(
    adata,
    tree_yaml: str,
    markers_yaml: str,
    root: str = 'root',
    min_marker_count: int = 2,
    verbose: bool = False,
    smoothing: bool = True,        # one-pass KNN score smoothing per split
    majority_vote: bool = True,    # optional single-pass label MV (excludes self)
    save_scores: bool = False,     # persist per-split cell-type scores as "<child>_score"
    min_score: float = 0.0,        # gate final labels below this raw score as "unknown"
    min_pruning_fc: float = 2.0    # minimum FC vs avg(other siblings) for +markers; set ≤1 to disable
):
    """
    TreeTag: hierarchical cell-type annotation using marker genes on a directed tree.

    What it does
    - Walks the tree top-down from `root` and assigns each cell to one child at every split.
    - Computes per-child scores from +/− marker averages and scales them robustly to [0,1] (1–99% quantiles).
    - (Optional) Prunes weak +markers by fold-change against the **average of the other siblings** at each split.
    - (Optional) Smooths scores with KNN neighbors and/or applies majority-vote on final labels.
    - (Optional) Gates low-confidence labels to "unknown" using a score threshold.

    Inputs
    - adata: AnnData with gene expression in .X (or .raw.X if present). If smoothing/majority_vote=True,
      adata must also contain a neighbor connectivities matrix in .obsp (e.g. 'connectivities').
    - tree_yaml, markers_yaml: YAML files consumed by init_tree(...), which attaches per-node marker lists.
    - root: **start node** for tagging (can be any existing vertex name in the tree).

    Outputs / side effects
    - Writes categorical labels to adata.obs["TreeTag"].
    - If save_scores=True, writes per-child score columns named "<child>_score" for splits that were scored.

    Performance & caching
    - Uses RAM-only caches keyed by this AnnData and the (sub)tree fingerprints; nothing is written to .uns or .h5ad.
    - Marker columns are pre-scaled by their nonzero max once per session; computations are sparse-friendly.
    """
    # Imports are scoped here so the function is self-contained.
    import time, hashlib, weakref
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix, csr_matrix, coo_matrix, issparse

    # Utility helpers used throughout TreeTag.
    def _row_normalize(M):
        """Row-stochastic normalization for CSR matrices; safe on zero rows."""
        M = M.tocsr()
        rs = np.asarray(M.sum(axis=1)).ravel()
        rs[rs == 0] = 1.0
        M = M.multiply((1.0 / rs)[:, None]).tocsr()
        return M

    def _quantile_01(A):
        """Per-column 1–99% scaling to [0,1]. Returns float32."""
        q1, q99 = np.quantile(A, [0.01, 0.99], axis=0)
        S = (A - q1) / (q99 - q1 + 1e-9)
        np.clip(S, 0, 1, out=S)
        return S.astype(np.float32, copy=False)

    def _log(msg: str):
        if verbose:
            print(msg)

    def _build_weights(used_cols, pos_lists, neg_lists, k=None):
        """Build W_pos/W_neg from marker index lists and chosen used_cols."""
        if k is None:
            k = len(pos_lists)
        local_map = {g: i for i, g in enumerate(used_cols)}
        data_p, rows_p, cols_p = [], [], []
        data_n, rows_n, cols_n = [], [], []
        for j, (pidx, nidx) in enumerate(zip(pos_lists, neg_lists)):
            if pidx:
                w = 1.0 / float(len(pidx))
                rows_p.extend(local_map[g] for g in pidx)
                cols_p.extend([j] * len(pidx))
                data_p.extend([w] * len(pidx))
            if nidx:
                w = 1.0 / float(len(nidx))
                rows_n.extend(local_map[g] for g in nidx)
                cols_n.extend([j] * len(nidx))
                data_n.extend([w] * len(nidx))
        n_local = len(used_cols)
        W_pos = coo_matrix((data_p, (rows_p, cols_p)), shape=(n_local, k)).tocsr()
        W_neg = coo_matrix((data_n, (rows_n, cols_n)), shape=(n_local, k)).tocsr() if data_n else None
        return W_pos, W_neg

    def _prune_debug(c, child_list, counts, means_used, union_pos, drop_cols, marker_union):
        """Verbose diagnostics for pruned markers at a split."""
        try:
            ct_str = ", ".join(f"{child_list[ii]}:{int(counts[ii])}" for ii in range(len(child_list)))
            _log(f"[Prune:DETAIL] {c}: cells per child -> " + ct_str)
            child_idx = child_list.index(c)
            for col_global in drop_cols:
                ui = union_pos[col_global]
                gene = marker_union[col_global]
                vals = [means_used[i, ui] for i in range(len(child_list))]
                sum_mu = np.nansum(vals)
                k = len(vals)
                other_mean = (sum_mu - vals[child_idx]) / max(1, k-1)
                child_val = vals[child_idx]
                EPS_dbg = 1e-9
                fc_dbg = (float(child_val) + EPS_dbg) / (float(other_mean) + EPS_dbg) if np.isfinite(child_val) else float('nan')
                vals_fmt = "; ".join(
                    f"{child_list[i]}={float(v):.4f}" if np.isfinite(v) else f"{child_list[i]}=nan"
                    for i, v in enumerate(vals)
                )
                other_str = f"{other_mean:.4f}" if np.isfinite(other_mean) else "nan"
                child_str = f"{child_val:.4f}" if np.isfinite(child_val) else "nan"
                _log(f"   - {gene}: {vals_fmt} | other_avg={other_str}, child={child_str}, FC(child/others)={fc_dbg:.3f}")
        except Exception as _e:
            _log(f"[Prune:DETAIL] (failed to print details: {_e})")

    # Session RAM cache: per-AnnData buckets keyed by id(adata), auto-cleaned when the AnnData is garbage-collected.
    global _TT_CACHE, _TT_CACHE_FINALIZERS
    try:
        _TT_CACHE
        _TT_CACHE_FINALIZERS
    except NameError:
        _TT_CACHE = {}
        _TT_CACHE_FINALIZERS = {}

    def _tt_cleanup(key):
        _TT_CACHE.pop(key, None)
        _TT_CACHE_FINALIZERS.pop(key, None)

    def _tt_get_bucket(adata):
        k = id(adata)
        bucket = _TT_CACHE.get(k)
        if bucket is None:
            bucket = {}
            _TT_CACHE[k] = bucket
            if k not in _TT_CACHE_FINALIZERS:
                _TT_CACHE_FINALIZERS[k] = weakref.finalize(adata, _tt_cleanup, k)
        return bucket

    # Ensure no legacy on-disk cache from older versions is present in .uns.
    if isinstance(getattr(adata, "uns", None), dict):
        adata.uns.pop("_treetag_cache", None)

    # Pruning threshold: set ≤1 to disable fold-change pruning of positive markers.
    thr = float(min_pruning_fc)
    _do_prune = np.isfinite(thr) and (thr > 1.0)
    _log(f"[TreeTag] pruning_min_fc={min_pruning_fc} → {'ON' if _do_prune else 'OFF'}")

    # Load the cell-type tree (init_tree attaches +/− markers to nodes). Note: root may be any node.
    use_raw = adata.raw is not None
    var_names = (adata.raw.var_names if use_raw else adata.var_names)

    G = init_tree(tree_yaml, markers_yaml, root=root, adata=adata)

    # Fast lookup: node name → graph index.
    vmap = {name: i for i, name in enumerate(G.vs["name"])}
    if root not in vmap:
        raise ValueError(f"Root '{root}' not found in tree.")

    # === Subtree support ===
    root_idx = G.vs.find(name=root).index
    subtree_idxs = G.subcomponent(root_idx, mode="OUT")
    subtree_vertices = [G.vs[i] for i in subtree_idxs]

    # Union of all +/− markers present anywhere in the SUBTREE only.
    marker_union = sorted({g for v in subtree_vertices for g in (v["pos_markers"] + v["neg_markers"])})

    # Retrieve (or create) this AnnData’s RAM cache bucket.
    cache = _tt_get_bucket(adata)

    # Hash fingerprints used to validate cached artifacts for this adata/(sub)tree combo.
    var_md5 = hashlib.md5(",".join(map(str, var_names)).encode("utf-8")).hexdigest()

    def _tree_fingerprint():
        parts = []
        for v in subtree_vertices:
            parts.append(v["name"])
            parts.append("|".join(v["pos_markers"]))
            parts.append("|".join(v["neg_markers"]))
        parts.append(f"ROOT::{root}")
        return hashlib.md5("\n".join(parts).encode("utf-8")).hexdigest()

    tree_fp = _tree_fingerprint()

    # Pre-scaled marker expression matrix (CSC), cached in RAM.
    def _get_scaled_marker_csc():
        """
        Return CSC (float32) of subtree marker union, scaled per-column by nonzero max.
        Cached in RAM via tokens; never persisted.
        """
        token = ("X_csc_v3",
                 int(adata.n_obs),
                 bool(use_raw),
                 tuple(marker_union),
                 var_md5,
                 tree_fp)
        entry = cache.get("X_csc_entry")
        if entry and entry.get("token") == token:
            return entry["X_csc"]

        # Extract the marker columns from X (or raw.X).
        col_map = {g: i for i, g in enumerate(var_names)}
        cols = [col_map[g] for g in marker_union if g in col_map]
        if not cols:
                    col_map = {g: i for i, g in enumerate(var_names)}
        cols = [col_map[g] for g in marker_union if g in col_map]
        if not cols:
            _log("[TreeTag] No marker overlap with current gene IDs; trying auto-convert to symbols...")
            try:
                info = convert(adata)  # uses var['feature_name']/['gene_symbols']/['SYMBOL'] if present
                # refresh gene list after conversion
                var_names = (adata.raw.var_names if use_raw else adata.var_names)
                col_map = {g: i for i, g in enumerate(var_names)}
                cols = [col_map[g] for g in marker_union if g in col_map]
                if not cols:
                    raise RuntimeError("After auto-conversion, still zero overlap between markers and genes.")
                _log(f"[TreeTag] Auto-converted genes (used={info.get('used')}, changed={info.get('changed')}).")
            except Exception as e:
                raise RuntimeError(
                    "No marker genes overlap your dataset and auto-conversion failed.\n"
                    "If your IDs are Ensembl, install gprofiler-official and convert first:\n"
                    "    pip install gprofiler-official\n"
                    "Then run: from treetag import convert; convert(adata)\n"
                    f"Details: {e}"
                )

        Xsrc = (adata.raw.X if use_raw else adata.X)
        t0 = time.perf_counter()
        X_mark = Xsrc[:, cols]
        X_csc = X_mark if (issparse(X_mark) and getattr(X_mark, "format", "") == "csc") else csc_matrix(X_mark)
        if X_csc.dtype != np.float32:
            X_csc = X_csc.astype(np.float32)

        # Scale each marker column by its nonzero max (sparse-friendly, done in place).
        indptr, data = X_csc.indptr, X_csc.data
        m = X_csc.shape[1]
        mx = np.zeros(m, dtype=np.float32)
        for j in range(m):
            seg = data[indptr[j]:indptr[j+1]]
            if seg.size:
                mx[j] = seg.max()
        nz = mx > 1e-12
        if nz.any():
            inv = np.zeros_like(mx)
            inv[nz] = 1.0 / mx[nz]
            for j in np.where(nz)[0]:
                data[indptr[j]:indptr[j+1]] *= inv[j]

        cache["X_csc_entry"] = {"token": token, "X_csc": X_csc}
        msg = "**pre-scaled**" if nz.any() else "pre-scaled (no scaling needed)"
        _log(f"[TreeTag] prepared {msg} marker CSC in {time.perf_counter()-t0:.2f}s (cells={X_csc.shape[0]}, genes={X_csc.shape[1]})")
        return X_csc

    X_csc = _get_scaled_marker_csc()
    g2col = {g: i for i, g in enumerate(marker_union)}

    def _children(node):
        # igraph one-liner: direct children of `node` (works for any subtree)
        return [G.vs[i]["name"] for i in G.successors(G.vs.find(name=node))]

    # Neighbor smoothing matrix: row-normalized KNN connectivities, cached in RAM.
    C_full = None
    if smoothing or majority_vote:
        def _neighbors_matrix():
            """Get normalized neighbor connectivities with caching."""
            conn_key = adata.uns.get("neighbors", {}).get("connectivities_key")
            if conn_key and conn_key in adata.obsp:
                C_raw = adata.obsp[conn_key]; src = conn_key
            elif "connectivities" in adata.obsp:
                C_raw = adata.obsp["connectivities"]; src = "connectivities"
            else:
                # Build neighbors via Scanpy if missing
                try:
                    import scanpy as sc
                except ImportError as e:
                    raise RuntimeError(
                        "smoothing/majority_vote requires a neighbor graph, and scanpy is missing. "
                        "Install it or set smoothing=False and majority_vote=False."
                    ) from e
                # Use whatever you prefer; X_pca if present, else .X
                use_rep = "X_pca" if "X_pca" in adata.obsm_keys() else None
                sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep)
                src = adata.uns.get("neighbors", {}).get("connectivities_key", "connectivities")
                C_raw = adata.obsp[src]
            nb_token = ("neighbors_v1", src, C_raw.shape, int(getattr(C_raw, "nnz", C_raw.size)))
            nb_entry = cache.get("neighbors_entry")
            if nb_entry and nb_entry.get("token") == nb_token:
                return nb_entry["C_full"]
            if C_raw.shape[0] != adata.n_obs:
                raise ValueError("neighbor connectivities shape mismatch with adata.n_obs")
            C_full_local = _row_normalize(csr_matrix(C_raw))
            cache["neighbors_entry"] = {"token": nb_token, "C_full": C_full_local}
            return C_full_local
        C_full = _neighbors_matrix()

    # Per-split cache: columns and weight matrices used for each tree node.
    node_cache_key = ("nodecache_v1", tree_fp, int(min_marker_count), tuple(marker_union))
    node_cache = cache.get("node_cache")
    if not (isinstance(node_cache, dict) and node_cache.get("_key") == node_cache_key):
        node_cache = {"_key": node_cache_key}
        cache["node_cache"] = node_cache

    def _prepare_node(node, kids):
        if node in node_cache:
            return node_cache[node]
        used_cols, seen = [], set()
        child_list, child_score_names, pos_lists, neg_lists = [], [], [], []
        for c in kids:
            v = G.vs[vmap[c]]
            pidx = [g2col[g] for g in v["pos_markers"] if g in g2col]
            nidx = [g2col[g] for g in v["neg_markers"] if g in g2col]
            if len(pidx) < max(1, min_marker_count):
                continue
            child_list.append(c)
            child_score_names.append(f"{c}_score")
            pos_lists.append(pidx)
            neg_lists.append(nidx)
            for gidx in pidx:
                if gidx not in seen: seen.add(gidx); used_cols.append(gidx)
            for gidx in nidx:
                if gidx not in seen: seen.add(gidx); used_cols.append(gidx)

        k = len(child_list)
        if k == 0:
            node_cache[node] = dict(used_cols=[], child_list=[], child_score_names=[], W_pos=None, W_neg=None)
            return node_cache[node]

        W_pos, W_neg = _build_weights(used_cols, pos_lists, neg_lists, k=k)

        node_cache[node] = dict(used_cols=used_cols, child_list=child_list,
                                child_score_names=child_score_names, W_pos=W_pos, W_neg=W_neg)
        return node_cache[node]

    # ---- NEW PRUNING: +markers must beat avg(other siblings) by FC >= thr ----
    def _prune_pos_markers_min_fc(row_idx, child_list, child_assignments):
        """
        Returns:
            kept_pos: dict(child -> list[int global_col_idx kept])
            pruned_log: dict(child -> list[str gene names pruned])
        """
        if (not _do_prune) or len(child_list) <= 1 or row_idx.size == 0:
            kept = {c: [g2col[g] for g in G.vs[vmap[c]]["pos_markers"] if g in g2col] for c in child_list}
            return kept, {c: [] for c in child_list}

        pos_by_child = {c: [g2col[g] for g in G.vs[vmap[c]]["pos_markers"] if g in g2col] for c in child_list}
        union_cols = sorted({g for lst in pos_by_child.values() for g in lst})
        if not union_cols:
            _log("[Prune] No positive markers present at this split; skipping pruning.")
            kept = {c: [] for c in child_list}
            return kept, {c: [] for c in child_list}

        # Group means on PRE-SCALED marker matrix
        k = len(child_list)
        X_rows = X_csc[row_idx, :]
        Xg = X_rows[:, union_cols]
        counts = np.array([np.count_nonzero(child_assignments == j) for j in range(k)], dtype=int)
        means_used = np.zeros((k, len(union_cols)), dtype=np.float64)
        for j in range(k):
            idx = np.where(child_assignments == j)[0]
            if idx.size == 0:
                continue
            Xj = Xg[idx, :]
            sums = np.asarray(Xj.sum(axis=0)).ravel()
            means_used[j, :] = sums / float(idx.size)
        means_used = np.where(np.isfinite(means_used), means_used, 0.0)

        # Average-of-others per child per gene
        sum_mu = np.nansum(means_used, axis=0)             # (g,)
        union_pos = {col: i for i, col in enumerate(union_cols)}
        kept_pos, pruned_log = {}, {}
        hard_min = int(max(2, min_marker_count))  # enforce at least 2

        EPS = 1e-9
        for c_idx, c in enumerate(child_list):
            cols_c = pos_by_child[c]
            if counts[c_idx] == 0 or not cols_c:
                kept_pos[c] = cols_c[:]
                pruned_log[c] = []
                why = "0 cells assigned" if counts[c_idx] == 0 else "no candidate +markers present"
                _log(f"[Prune] {c}: skipped pruning ({why}). Kept {len(cols_c)} markers.")
                continue

            pos_idx_in_union = [union_pos[col] for col in cols_c]
            child_mean = means_used[c_idx, pos_idx_in_union]                    # (m,)
            others_mean = (sum_mu[pos_idx_in_union] - child_mean) / max(1, k-1) # (m,)
            fc = (child_mean + EPS) / (others_mean + EPS)                       # child / avg(others)
            keep_mask = fc >= thr

            kept_cols = [col for col, keep in zip(cols_c, keep_mask) if keep]
            drop_cols = [col for col, keep in zip(cols_c, keep_mask) if not keep]
            kept_pos[c] = kept_cols
            pruned_log[c] = [marker_union[col] for col in drop_cols]

            _log(f"[Prune] {c}: kept {len(kept_cols)}/{len(cols_c)} (+markers, FC vs avg(others) ≥ {min_pruning_fc:g}).")
            if pruned_log[c]:
                _log("         pruned → " + ", ".join(pruned_log[c]))

            if len(kept_cols) < hard_min:
                msg = (f"[Prune:STOP] {c} has {len(kept_cols)} positive markers after pruning "
                       f"(min required = {hard_min}). Lower min_pruning_fc or review markers.")
                _log(msg)
                _prune_debug(c, child_list, counts, means_used, union_pos, drop_cols, marker_union)
                raise ValueError(msg)

        return kept_pos, pruned_log

    # Score children at a split (optionally prune and/or smooth).
    def _score_children_rows(row_idx, node, kids):
        cache_node = _prepare_node(node, kids)
        used_cols = cache_node["used_cols"]
        child_list = cache_node["child_list"]
        child_score_names = cache_node["child_score_names"]
        W_pos, W_neg = cache_node["W_pos"], cache_node["W_neg"]

        if not child_list:
            return [], [], np.zeros((row_idx.size, 0), dtype=np.float32), 0

        # First-pass scoring (no pruning)
        X_used = X_csc[:, used_cols][row_idx, :].tocsr()
        P = (X_used @ W_pos).toarray()
        N = (X_used @ W_neg).toarray() if W_neg is not None else None
        Diff = P if N is None else (P - N)
        np.clip(Diff, 0, None, out=Diff)

        # 1–99% split-wise quantile normalization → [0,1]
        S = _quantile_01(Diff)

        # Optional pruning based on this first pass
        if _do_prune and len(child_list) > 1:
            best = S.argmax(axis=1)
            kept_pos, _ = _prune_pos_markers_min_fc(row_idx, child_list, best)

            pos_lists_pruned = [kept_pos[c] for c in child_list]
            neg_lists = [[g2col[g] for g in G.vs[vmap[c]]["neg_markers"] if g in g2col] for c in child_list]

            seen, used_cols2 = set(), []
            for lst in pos_lists_pruned + neg_lists:
                for gidx in lst:
                    if gidx not in seen:
                        seen.add(gidx); used_cols2.append(gidx)

            if not used_cols2:
                return child_list, child_score_names, np.zeros((row_idx.size, 0), dtype=np.float32), 0

            W_pos2, W_neg2 = _build_weights(used_cols2, pos_lists_pruned, neg_lists, k=len(child_list))

            # re-score with pruned weights (X_csc already pre-scaled)
            X_used2 = X_csc[:, used_cols2][row_idx, :].tocsr()
            P2 = (X_used2 @ W_pos2).toarray()
            N2 = (X_used2 @ W_neg2).toarray() if W_neg2 is not None else None
            Diff2 = P2 if N2 is None else (P2 - N2)
            np.clip(Diff2, 0, None, out=Diff2)

            S = _quantile_01(Diff2)   # replace S
            used_cols = used_cols2                  # for logging

        # Optional smoothing (applies to either initial or pruned S)
        if smoothing and C_full is not None and S.shape[1] > 0:
            C_sub = C_full[row_idx, :][:, row_idx].tocsr()
            C_sub.setdiag(0.0)  # exclude self
            C_sub = _row_normalize(C_sub)
            S = C_sub @ S

        return child_list, child_score_names, S, len(used_cols)

    # Top-down traversal from the (sub)root: assign cells recursively.
    labels = pd.Series(index=adata.obs_names, dtype=object)
    n_obs = adata.n_obs
    winner_score = np.full(n_obs, np.nan, dtype=np.float32)
    created_cols = set()

    def _save_split_scores(score_names, row_idx, S):
      """Persist per-child scores, even if no cells were assigned there."""
      if not save_scores or not score_names:
          return
      for j, col in enumerate(score_names):
          if col not in created_cols:
              # Preallocate column for *all* cells with NaN (instead of zeros)
              adata.obs[col] = np.full(n_obs, np.nan, dtype=np.float32)
              created_cols.add(col)
          # If S has this column, write scores for the relevant rows
          if j < S.shape[1]:
              adata.obs[col].values[row_idx] = S[:, j]
          else:
              # No scores computed for this child at this split → leave as NaN
              adata.obs[col].values[row_idx] = np.nan

    def _walk(node, row_idx):
        kids = _children(node)
        if not kids or row_idx.size == 0:
            labels.iloc[row_idx] = node
            return
        t0 = time.perf_counter()
        child_names, child_score_names, S, n_local = _score_children_rows(row_idx, node, kids)
        t1 = time.perf_counter()
        _log(f"[TreeTag] split @ '{node}' (cells={row_idx.size}, X_used={row_idx.size}x{n_local}, score={t1 - t0:.2f}s, children_scored={len(child_names)})")
        if S.shape[1] == 0:
            labels.iloc[row_idx] = node
            return
        _save_split_scores(child_score_names, row_idx, S)
        best = S.argmax(axis=1)
        for j, child in enumerate(child_names):
            mask = (best == j)
            if mask.any():
                winner_score[row_idx[mask]] = S[mask, j]
                _walk(child, row_idx[mask])

    t_start = time.perf_counter()
    _walk(root, np.arange(adata.n_obs, dtype=int))
    _log(f"[TreeTag] total {time.perf_counter() - t_start:.2f}s")

    # Write final labels to adata.obs["TreeTag"].
    labels = labels.fillna("unidentified")
    cats = pd.Index(pd.unique(labels.values))
    adata.obs["TreeTag"] = pd.Categorical(labels.values, categories=cats)

    # Optional majority vote: one-pass neighbor label smoothing.
    if majority_vote:
        if C_full is None:
            raise ValueError("majority_vote requested but neighbor connectivities were not prepared")
        label_to_idx = {lab: i for i, lab in enumerate(cats)}
        idx = np.fromiter((label_to_idx.get(lab, -1) for lab in adata.obs["TreeTag"].astype(str).values), dtype=int)
        rows = np.where(idx >= 0)[0]
        cols = idx[idx >= 0]
        data = np.ones(rows.size, dtype=np.float32)
        OH = csr_matrix((data, (rows, cols)), shape=(adata.n_obs, len(cats)))

        C_mv = C_full.tocsr(copy=True)
        C_mv.setdiag(0.0)  # exclude self
        C_mv = _row_normalize(C_mv)

        votes = C_mv @ OH
        best_idx = votes.argmax(axis=1).A1
        new_labels = np.array(cats)[best_idx]
        adata.obs["TreeTag"] = pd.Categorical(new_labels, categories=cats)

    # Optional confidence gate: relabel low-score cells as "unknown".
    if (min_score is not None) and (float(min_score) > 0.0):
        final_labels = adata.obs["TreeTag"].astype(object).to_numpy()
        n = final_labels.size
        chosen = np.full(n, np.nan, dtype=np.float32)

        # use persisted per-label scores where available
        for lab in pd.unique(final_labels):
            if lab is None:
                continue
            col = f"{lab}_score"
            if col in adata.obs:
                vals = adata.obs[col].to_numpy(dtype=np.float32, copy=False)
                m = (final_labels == lab)
                chosen[m] = vals[m]

        # fallback to raw winner_score where needed
        missing = ~np.isfinite(chosen)
        if missing.any():
            tmp = winner_score.copy()
            tmp[~np.isfinite(tmp)] = -np.inf
            chosen[missing] = tmp[missing]

        # apply unknown threshold
        unknown_mask = (chosen < float(min_score))
        if np.any(unknown_mask):
            new_vals = final_labels.copy()
            new_vals[unknown_mask] = "unknown"
            cats = list(adata.obs["TreeTag"].cat.categories)
            if "unknown" not in cats:
                cats.append("unknown")
            adata.obs["TreeTag"] = pd.Categorical(new_vals, categories=cats)
