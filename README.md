# TreeTag

TreeTag is a lightweight Python package that automatically annotates single-cell RNA-seq data. It reads two editable YAML files: one lays out the hierarchy of cell types, and the other lists positive and negative marker genes. TreeTag promotes quick, interactive adjustment of marker sets and ontologies by keeping marker rules human-readable, and performing near-instant reannotation. Marker pruning avoids misleading assignments from dataset- or batch-specific markers, while smoothing helps overcome inherent scRNA-seq sparsity by integrating consistent signals from a PCA-driven neighborhood embedding.

---

## Features

- Integrates smoothly with AnnData and scanpy.

- Reads human‑editable YAMLs for the ontology and for positive/negative markers and builds the ontology as a graph (via igraph) for hierarchical traversal **(init_tree)**.

- Visualizes the ontology to inspect and validate structure **(plot_tree)**.

- Pre‑scales marker columns (sparse‑friendly), cache, and run lean matrix operations for fast scoring. **(TreeTag)**

- Computes hierarchical marker‑based scores top‑down; optionally applying KNN smoothing and majority vote using a PCA‑driven neighborhood embedding. **(TreeTag toggles)**

- Assigns cell‑type tags **(TreeTag in AnnData object)**.

- Prunes unreliable markers when they fail to separate the intended type  **(TreeTag toggles)**.

- Exposes per‑cell scores for manual inspection within AnnData/scanpy. **(*_score in AnnData object)**.

- Detects likely doublets after scoring, using per‑node scores to flag candidates for review/removal **(find_doublets)**.

---

## Installation
From PyPI (recommended)
```bash
pip install treetag
```
Upgrade
```bash
pip install --upgrade treetag
```
Verify intallation
```bash
python -c "import treetag, sys; print('TreeTag', treetag.__version__)"
```
---
## Quickstart

```python
import scanpy as sc
from treetag import TreeTag

# 1) Load data (PBMCs if you want it to work wuth the example YAML files)
adata = sc.read_h5ad("my_data.h5ad")

# 2) Prepare neighbors (required for smoothing / majority_vote)
sc.pp.pca(adata)
sc.pp.neighbors(adata, use_rep="X_pca")

# 3) Import example YAML files
data_dir = files("treetag.data")           # package submodule with YAMLs
tree_yaml = data_dir / "PBMC.yaml"
markers_yaml = data_dir / "PBMC_markers.yaml"

# 4) Run TreeTag
TreeTag(
    adata,
    tree_yaml="PBMC.yaml",     # cell-type hierarchy
    markers_yaml="PBMC_markers.yaml",   # positive/negative markers
    root="PBMC",                   # any node in your ontology
    save_scores=True                # optional: write per-node scores
)

# 5) Inspect results
print(adata.obs["TreeTag"].value_counts())
sc.pl.umap(adata, color="TreeTag")
```

## YAML File Formats

#### Ontology YAML

```yaml
root:
  T_NK:
    CD4_T:
      Treg:
      Th:
      _Tfh:
    CD8_T:
  B:
    Naive_B:
    Memory_B:
  Myeloid:
    Mono:
    DC:
```

**`!` note:** Any key starting with "_" is treated as **disabled**; the cell-type and its entire subtree are skipped.

#### Markers YAML

```yaml
T_NK: [CD2, IL32, CD7, CD247, CD3E, LCK, IFITM1, GIMAP7, -MS4A1]
CD4_T: [CD4, TRAT1, ICOS, GPR183, CD40LG, IL6ST, -CD8A, -CD8B]
Treg: [FOXP3, RTKN2, IL2RA, IKZF2, CTLA4, TNFRSF18, TIGIT, -CD40LG]
```
**`!` note:** At least 2 positive markers are needed per cell type. Negative markers start with "-" and are not obligatory. 

## Function reference

## `TreeTag`

**What it does:** Hierarchical cell‑type tagging using positive/negative markers.

**Signature:**

```python
TreeTag(
    adata,
    tree_yaml: str,
    markers_yaml: str,
    root: str,
    min_marker_count: int = 2,
    verbose: bool = False,
    smoothing: bool = True,
    majority_vote: bool = True,
    save_scores: bool = False,
    min_score: float = 0.0,
    min_pruning_fc: float = 1.5,
) -> None
```

**Writes:** `adata.obs["TreeTag"]`; if `save_scores=True`, also `<node>_score` columns.

**Requires (if enabled):** neighbors in `adata.obsp` for `smoothing`/`majority_vote`.

**Common errors (and fixes):**

* *No neighbor graph:* run `sc.pp.neighbors(adata, use_rep="X_pca")` **or** set `smoothing=False, majority_vote=False`.
* *No subtree markers found:* check gene naming (symbols vs Ensembl), root, and `.raw` usage.
* *Neighbor shape mismatch:* rebuild neighbors **after** any cell filtering.

---

### `init_ontology`

**What it does:** Loads ontology + markers, builds the graph, normalizes marker dicts. If adata is provided, omits missing markers.

**Signature:**

```python
G = init_ontology(tree_yaml: str, markers_yaml: str, root: str, adata=None,)
```

**Returns:**

* `G`: graph of the ontology (node names in `G.vs["name"]`), with poaitive and negative marker attributes per node.

---

### `markers`

**What it does:** Returns marker genes for a node (optionally filtered to genes present in `adata`).

**Signature:**

```python
markers(
    cell_type: str,
    sign: str = "pos",            # "pos" or "neg"
    markers_yaml: str = "markers.yaml",
    tree_yaml: str = "ontology.yaml",
    adata=None,                    # optional filter to adata.var_names/raw.var_names
) -> list[str]
```

---

### `subscores`

**What it does:** Lists existing `<node>_score` columns under a root (useful after `TreeTag(save_scores=True)`).

**Signature:**

```python
subscores(
    root_cell: str,
    adata,
    markers_yaml: str,
    tree_yaml: str,
) -> list[str]
```

---

### `find_doublets`

**What it does:** Flags likely doublets **after scoring** using per‑node score patterns (e.g., strong scores for incompatible lineages).

**Signature (minimal):**

```python
find_doublets(
    adata,
    threshold: float = 0.25,   # heuristic overlap metric; implementation‑specific
    write: bool = True,
    key: str = "doublet_like",
) -> "pd.Series[bool] | np.ndarray[bool]"
```

**Writes (if `write=True`):** `adata.obs["doublet_like"]` boolean mask.

---

### `plot_tree`

**What it does:** Renders the ontology tree (optionally overlaying counts/assignments).

**Signature (typical):**

```python
plot_tree(
    tree_yaml: str | None = None,
    markers_yaml: str | None = None,
    root: str | None = None,
    G=None,                      # alternatively pass a prebuilt graph
    adata=None,                  # optional: color by counts/labels
    ax=None,
    layout: str = "rt",         # e.g., top‑down
) -> "matplotlib.axes.Axes"
```

