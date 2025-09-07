# tests/test_markers_neg_both.py
import treetag as tt


def test_markers_neg_and_both(tiny_adata_csr, toy_yaml_pair):
    _, mark_y = toy_yaml_pair
    both = tt.markers(
        cell_type="B", sign="both", markers_yaml=mark_y, adata=tiny_adata_csr
    )
    pos = tt.markers(
        cell_type="B", sign="pos", markers_yaml=mark_y, adata=tiny_adata_csr
    )
    neg = tt.markers(
        cell_type="B", sign="neg", markers_yaml=mark_y, adata=tiny_adata_csr
    )

    assert isinstance(both, dict) and set(both.keys()) == {"pos", "neg"}
    assert set(pos) == {"MS4A1", "IGHM"}
    assert set(neg) == {"CD3D"}
    assert set(both["pos"]) == set(pos)
    assert set(both["neg"]) == set(neg)
