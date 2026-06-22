"""
Watch change-detection + live.html writing.
Run with pytest, or `python tests/test_watch.py`.
"""
import os
import tempfile

from henrietta_slit_guiding import watch


def test_decide_reason():
    fp, ev, cf = (10, 100, 1), 1000.0, 500.0
    assert watch.decide_reason(fp, ev, cf, None, None, None) == "startup"
    assert watch.decide_reason(fp, ev, cf, fp, ev, cf) is None
    assert "frame" in watch.decide_reason((11, 101, 1), ev, cf, fp, ev, cf)
    assert watch.decide_reason(fp, 2000.0, cf, fp, ev, cf) == "motion_events.txt edit"
    assert watch.decide_reason(fp, ev, 999.0, fp, ev, cf) == "config.txt edit"
    r = watch.decide_reason((11, 101, 1), 2000.0, cf, fp, ev, cf)
    assert "frame" in r and "events" in r
    # newest frame still being written (size grows, same max) -> redraw
    assert watch.decide_reason((10, 100, 99), ev, cf, fp, ev, cf) is not None


def test_ensure_live_html():
    d = tempfile.mkdtemp()
    watch.ensure_live_html(d)
    html = open(os.path.join(d, "live.html")).read()
    assert "dx_dy_vs_frame.png" in html and "box_overlay.png" in html


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
