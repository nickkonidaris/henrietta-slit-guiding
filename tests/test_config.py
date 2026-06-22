"""
Config core: night-date rule, run-dir naming, parsing + fallbacks, bundled
data, seeded template round-trip.  Run with pytest, or `python tests/test_config.py`.
"""
import datetime
import os
import tempfile

from henrietta_slit_guiding import config as c


def _tmp(text):
    fd, path = tempfile.mkstemp(suffix=".txt")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


def test_night_date_noon_rule():
    assert c.night_date(datetime.datetime(2026, 6, 22, 1, 30)) == "2026_06_21"
    assert c.night_date(datetime.datetime(2026, 6, 22, 11, 59)) == "2026_06_21"
    assert c.night_date(datetime.datetime(2026, 6, 21, 12, 0)) == "2026_06_21"
    assert c.night_date(datetime.datetime(2026, 6, 21, 23, 30)) == "2026_06_21"


def test_run_dir_name():
    n = c.run_dir_name("HD136352", "nu01Lup", datetime.datetime(2026, 6, 21, 20, 0))
    assert n == "2026_06_21-HD136352-nu01Lup"
    # spaces sanitized
    assert " " not in c.run_dir_name("nu 01", "x", datetime.datetime(2026, 6, 21, 20, 0))


def test_bundled_data():
    assert os.path.exists(c.bundled_bpm_path())
    assert "source" in c.template_text()


def test_parse_good_and_roles():
    p = _tmp("source /data/foo\nbpm bundled\nstart_frame 83\nobject 136352\n"
             "filter R-J\npa 9.34\npixscale 0.5\n"
             "box StarA target 100 20 600 800\nbox StarB comp 1900 20 600 800\n")
    cfg = c.load_config(p)
    assert cfg["source"] == "/data/foo"
    assert cfg["start_frame"] == 83 and cfg["object"] == "136352"
    assert cfg["filter"] == "R-J" and abs(cfg["pa"] - 9.34) < 1e-9
    assert cfg["bpm_path"] == c.bundled_bpm_path()
    t, comp = c.target_comp(cfg["boxes"])
    assert t == "StarA" and comp == "StarB"
    b = cfg["boxes"]["StarA"]
    assert b["x_lo"] == 80 and b["x_hi"] == 121
    os.unlink(p)


def test_bpm_relative_resolved_to_config_dir():
    p = _tmp("bpm sub/mask.fits\nbox A target 1 1 1 2\nbox B comp 2 1 1 2\n")
    cfg = c.load_config(p)
    assert cfg["bpm_path"] == os.path.join(os.path.dirname(os.path.abspath(p)), "sub/mask.fits")
    os.unlink(p)


def test_fallbacks():
    def is_default(cfg):
        d = c.DEFAULT_BOXES
        return set(cfg["boxes"]) == set(d)
    assert is_default(c.load_config("/no/such/file.txt"))
    for txt in ["box A target 1 1\n",                                  # too few
                "box A target 1 1 1 2\nbox B target 2 1 1 2\n",        # two targets
                "box A target 1 -5 1 2\nbox B comp 2 1 1 2\n"]:        # bad geom
        p = _tmp(txt)
        assert is_default(c.load_config(p))
        os.unlink(p)


def test_seeded_config_roundtrips():
    p = _tmp(c.seeded_config_text("HD1", "nuX"))
    cfg = c.load_config(p)
    t, comp = c.target_comp(cfg["boxes"])
    assert t == "HD1" and comp == "nuX"
    assert cfg["object"] == "HD1"
    os.unlink(p)


def test_update_config_boxes_update_preserves_other_lines():
    p = _tmp("source /data/foo\n# a comment\n"
             "box StarA target 100 20 600 800\nbox StarB comp 1900 20 600 800\n")
    rep = dict(c.update_config_boxes(p, {
        "StarA": dict(role="target", x_center=111, halfwidth=20, y_lo=600, y_hi=800),
        "StarB": dict(role="comp", x_center=1950, halfwidth=20, y_lo=600, y_hi=800),
    }))
    assert rep == {"StarA": "updated", "StarB": "updated"}
    assert os.path.exists(p + ".bak")
    cfg = c.load_config(p)
    assert cfg["source"] == "/data/foo"                 # untouched line preserved
    assert cfg["boxes"]["StarA"]["x_lo"] == 111 - 20    # updated center
    assert cfg["boxes"]["StarB"]["x_lo"] == 1950 - 20
    os.unlink(p); os.unlink(p + ".bak")


def test_update_config_boxes_appends_missing():
    p = _tmp("source /data\nbox StarA target 100 20 600 800\n")   # missing comp
    rep = dict(c.update_config_boxes(p, {
        "StarB": dict(role="comp", x_center=1950, halfwidth=20, y_lo=600, y_hi=800)}))
    assert rep == {"StarB": "added"}
    cfg = c.load_config(p)                               # now valid: 1 target + 1 comp
    t, comp = c.target_comp(cfg["boxes"])
    assert t == "StarA" and comp == "StarB"
    os.unlink(p); os.unlink(p + ".bak")


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
