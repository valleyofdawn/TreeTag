import scanpy as sc
from treetag import TreeTag, subscores, find_doublets

adata = sc.read_h5ad(r"data.h5ad")  # TODO: put your file path

# First run: no neighbors required (keeps it simple)
TreeTag(
    adata,
    tree_yaml=r"ontology.yaml",      # TODO
    markers_yaml=r"markers.yaml",    # TODO
    root="B_cells",
    smoothing=False,
    majority_vote=False,
    save_scores=True,
)

# Collect existing subtree scores into a CSV
df = subscores("B_cells", adata, tree_yaml=r"ontology.yaml", return_df=True)
df.to_csv("B_cells_scores.csv", index=False)

# Optional: quick doublet diagnostic from first split
find_doublets(adata, tree_yaml=r"ontology.yaml", markers_yaml=None, root="B_cells")

# Save updated AnnData (labels in .obs["TreeTag"], scores in <child>_score)
adata.write("adata_treetag.h5ad")
print("Done.")
