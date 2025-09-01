from typing import Optional, Union, Sequence, Tuple
import textwrap
import matplotlib.pyplot as plt
import igraph as ig
 
def plot_tree(
    G: "ig.Graph",
    root: Optional[str] = None,
    vertex_size: int = 50,
    vertex_label_size: int = 9,
    bbox: Tuple[int, int] = (1400, 800),   # pixels
    margin: int = 50,
    palette: Optional[Union[dict, Sequence[str]]] = None,
    wrap_width: int = 6,
    dpi: Optional[int] = None,
    ax: Optional[plt.Axes] = None,
    show: bool = True,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot the ontology tree with igraph.

    palette:
      - None                 -> all vertices gray
      - dict {node: color}   -> per-node colors (others default gray)
      - list/tuple           -> one color per vertex (len == number of vertices)

    Returns (fig, ax).
    """
    # --- labels ---
    names = G.vs["name"] if "name" in G.vs.attributes() else [str(i) for i in range(G.vcount())]
    labels = [textwrap.fill(n.replace("_", " "), width=wrap_width, break_long_words=True) for n in names]

    # --- layout (top-down Reingold–Tilford) ---
    if root is not None:
        try:
            ridx = G.vs.find(name=root).index
        except ValueError as e:
            raise ValueError(f"root '{root}' not found in graph node names") from e
        layout = G.layout_reingold_tilford(root=[ridx], mode="out")
    else:
        layout = G.layout_reingold_tilford(mode="out")
    coords = [(x, -y) for x, y in layout]  # flip Y for top-down

    # --- colors ---
    default_color = "#C0C0C0"
    if palette is None:
        vcols = [default_color] * len(names)
    elif isinstance(palette, dict):
        vcols = [palette.get(n, default_color) for n in names]
    elif isinstance(palette, (list, tuple)):
        if len(palette) != len(names):
            raise ValueError(f"palette length {len(palette)} != number of vertices {len(names)}")
        vcols = list(palette)
    else:
        raise TypeError("palette must be None, dict, or list/tuple of colors")

    # --- figure/axes ---
    if ax is None:
        # Figure size in inches
        use_dpi = dpi if dpi is not None else plt.rcParams.get("figure.dpi", 100)
        figsize = (bbox[0] / use_dpi, bbox[1] / use_dpi)
        fig, ax = plt.subplots(figsize=figsize, dpi=use_dpi)
    else:
        fig = ax.figure
    ax.set_axis_off()

    # --- draw ---
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
    if show:
        plt.show()
    return fig, ax
