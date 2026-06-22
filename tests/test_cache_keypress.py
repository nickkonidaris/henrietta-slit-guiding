"""
Fit-cache signature/round-trip/invalidation and keypress matrix math.
Run with pytest, or `python tests/test_cache_keypress.py`.
"""
import os
import tempfile

import numpy as np

from henrietta_slit_guiding import config as c
from henrietta_slit_guiding import montage, keypress

CFG = {"bpm_path": "/x/bpm.fits"}
BOXES = c.DEFAULT_BOXES


def test_signature_changes_on_template():
    a = montage.cache_signature(CFG, 83, BOXES)
    b = montage.cache_signature(CFG, 99, BOXES)
    assert a != b


def test_signature_changes_on_box():
    other = {k: dict(v) for k, v in BOXES.items()}
    first = next(iter(other))
    other[first]["x_lo"] -= 5
    assert montage.cache_signature(CFG, 83, BOXES) != montage.cache_signature(CFG, 83, other)


def test_cache_roundtrip_and_invalidation():
    sigA = montage.cache_signature(CFG, 83, BOXES)
    sigB = montage.cache_signature(CFG, 99, BOXES)
    frames = {"100": {"size": 123, "time_mjd": 60000.0,
                      "stars": {"TARGET": {"dx": 0.5, "dy": -0.2}}}}
    fd, path = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    montage.save_cache(sigA, frames, path)
    assert montage.load_cache(sigA, path) == frames
    assert montage.load_cache(sigB, path) == {}        # signature mismatch -> empty
    os.unlink(path)
    assert montage.load_cache(sigA, "/no/such.json") == {}


def test_keypress_matrix_invertible():
    M, _, _ = keypress.response_matrix(9.34, 1.0)
    assert abs(np.linalg.det(M)) > 1e-3


def test_keypress_deadband_and_presses():
    assert keypress.recommend(0.05, 0.05, 9.34)["within_deadband"] is True
    rec = keypress.recommend(1.3, -0.4, 9.34)
    assert rec["within_deadband"] is False
    assert isinstance(rec["presses"], list)
    # nulling the drift: applying the recommended exact shift should reduce |d|
    M, _, _ = keypress.response_matrix(9.34, 1.0)
    shift = M @ np.array([rec["n_up"], rec["n_left"]])
    resid = np.array([1.3, -0.4]) + shift
    assert np.hypot(*resid) < np.hypot(1.3, -0.4)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    fails = 0
    for fn in fns:
        try:
            fn()
            print(f"  [ok  ] {fn.__name__}")
        except Exception as e:  # noqa: BLE001
            fails += 1
            print(f"  [FAIL] {fn.__name__}: {e}")
    raise SystemExit(1 if fails else 0)
