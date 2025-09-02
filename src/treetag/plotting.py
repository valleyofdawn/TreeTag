import textwrap, colorsys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import igraph as ig
from ._init_tree import _init_tree as init_tree

def _hsv_palette(n, s=0.5, v=0.85):
    if n <= 0:
        return []
    return [mcolors.to_hex(colorsys.hsv_to_rgb(i / n, s, v)) for i in range(n)]

def plot_tree(
    tree_yaml,
    root="root",
    vertex_size=50,
    vertex_label_size=9,
    bbox=(1400, 800),   # pixels
    margin=50,
    palette=None,       # None -> auto HSV pastel; str -> matplotlib cmap name
    wrap_width=6,
):
    """
    Plot the ontology tree with igraph.

    palette:
      - None  -> auto HSV pastel (S=0.3, V=1)
      - str   -> matplotlib colormap name ("tab20", "turbo", "rainbow", ...)
    """
    G = init_tree(tree_yaml, root=root)
    names = G.vs["name"] if "name" in G.vs.attributes() else [str(i) for i in range(G.vcount())]
    labels = [textwrap.fill(n.replace("_", " "), width=wrap_width, break_long_words=True) for n in names]

    # layout (top-down)
    ridx = G.vs.find(name=root).index
    layout = G.layout_reingold_tilford(root=[ridx], mode="out")
    coords = [(x, -y) for x, y in layout]

    # colors
    n = len(names)
    if palette is None:
        vcols = _hsv_palette(n, s=0.3, v=1)
    else:
        cmap = plt.get_cmap(palette)
        vcols = [mcolors.to_hex(cmap(i / max(n - 1, 1))) for i in range(n)]

    # figure
    fig, ax = plt.subplots(figsize=(bbox[0] / 100, bbox[1] / 100))
    ax.set_axis_off()

    # draw
    ig.plot(
        G.as_undirected(),
        target=ax,
        layout=coords,
        vertex_label=labels,
        vertex_size=vertex_size,
        vertex_label_size=vertex_label_size,
        margin=margin,
        vertex_color=vcols,
        vertex_frame_color="dimgray",
        edge_color="gray",
    )
    plt.tight_layout()
    plt.show()
    return fig, ax
