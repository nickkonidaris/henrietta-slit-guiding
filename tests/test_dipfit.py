"""
Regression: the Y-fit window must bracket the companion's low absorption dip.

The companion's feature lands a few rows lower than the target's (isowavelength
tilt across the detector); with the window centered on box-row 30 ([15,45]) a
dip near row 23 must be fit at its true center, NOT railed to a window bound.
This guards against regressing GAUSS_GUESS_ROW back up (e.g. to 40 -> [25,55]),
which clipped the companion's dip.

Run with pytest, or `python tests/test_dipfit.py`.
"""
import numpy as np

from henrietta_slit_guiding import montage as m


def _spectrum_with_dip(center, n=182, amp=-3000.0, sigma=3.0):
    x = np.arange(n, dtype=float)
    return 100.0 + amp * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))


def test_companion_low_dip_not_railed():
    lo = m.GAUSS_GUESS_ROW - m.GAUSS_HALF_WINDOW
    hi = m.GAUSS_GUESS_ROW + m.GAUSS_HALF_WINDOW
    _, _, mu, _, _ = m.fit_dip(_spectrum_with_dip(23.0))
    assert abs(mu - 23.0) < 1.5, f"dip at 23 fit to mu={mu}"
    assert abs(mu - lo) > 1.0 and abs(mu - hi) > 1.0, f"mu={mu} railed to bound [{lo},{hi}]"


def test_target_dip_fit():
    _, _, mu, _, _ = m.fit_dip(_spectrum_with_dip(27.0))
    assert abs(mu - 27.0) < 1.5, f"dip at 27 fit to mu={mu}"


def test_window_reaches_low_rows():
    # the whole point of the fix: window covers the companion's row
    assert (m.GAUSS_GUESS_ROW - m.GAUSS_HALF_WINDOW) <= 20


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
