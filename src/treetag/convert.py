# src/treetag/convert.py
from __future__ import annotations
from typing import Sequence


def convert(
    adata,
    prefer_var_cols: Sequence[str] = ("feature_name", "gene_symbols", "SYMBOL"),
):
    """
    Minimal converter.

    - If a symbol column exists in `adata.var` (e.g. 'feature_name'), use it to set `adata.var_names`.
    - Stores originals in `adata.var['original_gene']` and, when shapes match, in `adata.raw.var['original_gene']`.
    - If no such column exists:
        * If var_names look like Ensembl IDs -> raise with an instruction to install gprofiler and convert.
        * Otherwise assume var_names are already symbols and do nothing.

    Returns: dict with {'changed': int, 'used': str} where used is the column name or 'none'.
    """
    import numpy as np
    import pandas as pd

    # ---- use existing symbol column, if present ----
    for col in prefer_var_cols:
        if col in adata.var.columns:
            vals = adata.var[col].astype(str)
            # require that it's mostly populated and not constant
            if vals.notna().mean() >= 0.5 and vals.nunique() > 10:
                orig = adata.var_names.astype(str)

                # rename source column first to avoid index/column clashes
                adata.var.rename(columns={col: f"{col}_orig"}, inplace=True)
                adata.var["original_gene"] = orig

                adata.var_names = vals
                adata.var_names_make_unique()

                # best-effort sync of .raw (only if same number of vars)
                if (
                    adata.raw is not None
                    and getattr(adata.raw, "n_vars", None) == adata.n_vars
                ):
                    try:
                        adata.raw.var["original_gene"] = adata.raw.var_names.astype(str)
                        adata.raw.var.index = adata.var_names
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
            "before running TreeTag. Example:\n"
            "    pip install gprofiler-official\n"
            "    # then run your own conversion step using g:Profiler\n"
        )

    # assume var_names are already symbols; do nothing
    return {"changed": 0, "used": "none"}
