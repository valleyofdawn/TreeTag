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

    def _disarm_index_column_collision(df: "pd.DataFrame") -> None:
        """If df.index.name matches a column with different values, rename that column."""
        idxn = df.index.name
        if idxn is not None and idxn in df.columns:
            # compare as strings for robustness
            col_vals = df[idxn].astype(str)
            idx_vals = pd.Index(df.index.astype(str))
            if not idx_vals.equals(pd.Index(col_vals)):
                df.rename(columns={idxn: f"{idxn}_orig"}, inplace=True)

    # Preempt collision on current .var
    _disarm_index_column_collision(adata.var)

    for col in prefer_var_cols:
        if col in adata.var.columns:
            vals = adata.var[col].astype(str)
            # require that it's mostly populated and not constant
            if vals.notna().mean() >= 0.5 and vals.nunique() > 10:
                orig = adata.var_names.astype(str)

                # rename source column first to avoid column re-use
                adata.var.rename(columns={col: f"{col}_orig"}, inplace=True)
                adata.var["original_gene"] = orig

                # set new var_names with name=None to avoid future collisions
                new_index = pd.Index(vals.values, name=None)
                adata.var_names = new_index
                adata.var.index.name = None  # belt-and-suspenders

                # best-effort sync of .raw (only if same number of vars)
                if (
                    adata.raw is not None
                    and getattr(adata.raw, "n_vars", None) == adata.n_vars
                ):
                    try:
                        rv = adata.raw.var
                        _disarm_index_column_collision(rv)  # guard raw.var too
                        rv["original_gene"] = adata.raw.var_names.astype(str)
                        rv.index = pd.Index(adata.var_names, name=None)
                    except Exception:
                        pass

                changed = int(np.sum(orig.values != adata.var_names.values))
                return {"changed": changed, "used": col}

    # ---- no column: decide whether to error or leave as-is ----
    import pandas as pd

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
