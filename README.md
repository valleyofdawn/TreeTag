# TreeTag

TreeTag is a lightweight Python package that automatically annotates single-cell RNA-seq data. It reads two editable YAML files: one lays out the hierarchy of cell types, and the other lists positive and negative marker genes. TreeTag promotes quick, interactive adjustment of marker sets and ontologies by keeping marker rules human-readable, and performing near-instant reannotation. Marker pruning avoids misleading assignments from dataset- or batch-specific markers, while smoothing helps overcome inherent scRNA-seq sparsity by integrating consistent signals from a PCA-driven neighborhood embedding.

---

## Features

- Automated scoring of marker gene sets
- Recursive splitting of heterogeneous populations
- Marker pruning and filtering
- Compatible with `scanpy` and `AnnData`
- Easily extendable and interpretable

---

## Installation

Coming soon as a pip-installable package. For now, clone the repository:

```bash
git clone https://github.com/yourusername/treetag.git
cd treetag
