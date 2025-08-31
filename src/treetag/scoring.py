from __future__ import annotations
import pandas as pd
from .tree import init_tree

def subscores(root_cell=root,
              adata,
              tree_yaml,
              only_leaves: bool = False):

    G  = init_tree(tree_yaml, markers_yaml, root=root_cell)

    root_v = G.vs.find(name=root_cell)

    # all descendants (OUT reachability) and drop the root itself
    desc_idxs = G.subcomponent(root_v.index, mode='OUT')
    desc_idxs = [i for i in desc_idxs if i != root_v.index]
    if only_leaves:
      desc_idxs = [i for i in desc_idxs if G.outdegree(i) == 0]

    nodes = [G.vs[i]['name'] for i in desc_idxs]

    existing, missing = [], []
    for n in nodes:
        col = f"{n}_score"
        (existing if col in adata.obs else missing).append(col)

    if missing:
        print(f"Missing score columns: {missing}")
    return existing
