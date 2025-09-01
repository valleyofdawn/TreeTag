# run_example.py
import anndata as ad
from treetag import TreeTag, subscores

adata = ad.read_h5ad(r"example.h5ad")

# first run: keep it simple (no neighbor graph needed)
TreeTag(
    adata,
    tree_yaml=r"data\PBMC_tree.yaml",     # <-- change if your file name is different
    markers_yaml=r"data\PBMC_markers.yaml",   # <-- change if your file name is different
    root="root",
    smoothing=False,
    majority_vote=False,
    save_scores=True,
)

# quick summary + outputs
print("TreeTag label counts:")
print(adata.obs["TreeTag"].value_counts(dropna=False))

# collect scores under the root subtree and save
df = subscores("root", adata, tree_yaml=r"ontology.yaml", return_df=True)
df.to_csv("example_scores.csv", index=False)

# save the labeled AnnData
adata.write("example_treetag.h5ad")
print("Wrote example_scores.csv and example_treetag.h5ad")
