# tests/test_smoke.py
def test_import():
    import treetag as tt

    assert hasattr(tt, "TreeTag")
