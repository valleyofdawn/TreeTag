__all__ = ["init_tree"]
__version__ = "0.1.0"

try:
    from .tree import init_tree
except Exception:
    # Allows packaging before you paste tree.py content
    def init_tree(*args, **kwargs):
        raise RuntimeError("init_tree not yet implemented")
