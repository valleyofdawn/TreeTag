# src/treetag/convert.py
from __future__ import annotations
from typing import Sequence


def convert(
    adata,
    prefer_var_cols: Sequence[str] = ("feature_name", "gene_symbols", "SYMBOL"),
):
    """
    Minimal converter with index/column collision guards.

    Returns: {'changed': int, 'used': str}
    """
    import numpy as np
    import pandas as pd

    def _disarm_index_column_collision(df: "pd.DataFrame") -> str | None:
        """
        If df.index.name matches a column with different values, rename that column.
        Returns the renamed column name if a rename occurred, else None.
        """
        idxn = df.index.name
        if idxn is not None and idxn in df.columns:
            col_vals = df[idxn].astype(str)
            idx_vals = pd.Index(df.index.astype(str))
            if not idx_vals.equals(pd.Index(col_vals)):
                new_col = f"{idxn}_orig"
                # avoid cascading renames
                suffix = 1
                while new_col in df.columns:
                    new_col = f"{idxn}_orig{suffix}"
                    suffix += 1
                df.rename(columns={idxn: new_col}, inplace=True)
                return new_col
        return None

    def _clear_index_name(df: "pd.DataFrame") -> None:
        # Clear index name to prevent future collisions when setting new indices
        if df.index.name is not None:
            df.index.name = None

    # Preempt collisions on current .var and .raw.var
    _disarm_index_column_collision(adata.var)
    _clear_index_name(adata.var)

    if adata.raw is not None:
        try:
            rv = adata.raw.var
            _disarm_index_column_collision(rv)
            _clear_index_name(rv)
        except Exception:
            pass

    # Helper: if a preferred column was just renamed due to collision, use the renamed version
    def _resolve_col(df: "pd.DataFrame", col: str) -> str | None:
        if col in df.columns:
            return col
        alt = f"{col}_orig"
        if alt in df.columns:
            return alt
        # Also handle numbered _orig variants from _disarm_index_column_collision
        candidates = [c for c in df.columns if c.startswith(f"{col}_orig")]
        return sorted(candidates, key=len)[0] if candidates else None

    for col in prefer_var_cols:
        src = _resolve_col(adata.var, col)
        if src is None:
            continue

        vals = adata.var[src].astype(str)
        # require that it's mostly populated and not constant
        if vals.notna().mean() >= 0.5 and vals.nunique() > 10:
            orig = adata.var_names.astype(str)

            # Preserve the source column by moving it out of the way if needed
            if src == col:
                keep_name = f"{col}_orig"
                if keep_name in adata.var.columns and keep_name != src:
                    # already preserved by disarm; do nothing
                    pass
                elif keep_name != src:
                    adata.var.rename(columns={src: keep_name}, inplace=True)
                    src = keep_name  # update handle

            # Track original symbols
            adata.var["original_gene"] = orig

            # Set new var_names with name=None to avoid collisions
            new_index = pd.Index(vals.values, name=None)
            adata.var_names = new_index  # may reorder .var index internally
            adata.var.index.name = None  # belt-and-suspenders

            # best-effort sync of .raw (only if same number of vars)
            if (
                adata.raw is not None
                and getattr(adata.raw, "n_vars", None) == adata.n_vars
            ):
                try:
                    rv = adata.raw.var
                    _disarm_index_column_collision(rv)
                    rv["original_gene"] = adata.raw.var_names.astype(str)
                    rv.index = pd.Index(adata.var_names, name=None)
                except Exception:
                    pass

            changed = int(np.sum(orig.values != adata.var_names.values))
            return {"changed": changed, "used": col}

    # ---- no column: decide whether to error or leave as-is ----
    vn = pd.Index(adata.var_names.astype(str))
    looks_ensembl = vn.str.startswith(("ENSG", "ENSMUSG", "ENSDARG")).mean() > 0.5

    if looks_ensembl:
        raise RuntimeError(
            "No symbol column found in `adata.var`, and gene IDs look like Ensembl. "
            "Please install gprofiler-official and perform ID conversion to gene symbols "
            "before running TreeTag."
        )

    # assume var_names are already symbols; do nothing
    return {"changed": 0, "used": "none"}
