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



