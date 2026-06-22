"""
The live loop: re-run the analysis whenever a new frame lands, or you edit
motion_events.txt or config.txt; write box_overlay.png at start + on config
change; serve an auto-refreshing live.html; and print a compact status
readout (drift + recommended keys) each redraw.
"""
from __future__ import annotations

import glob
import os
import sys
import time

import numpy as np

from . import config as cfgmod
from . import keypress
from . import montage

LIVE_HTML_NAME = "live.html"

LIVE_HTML_DOC = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>slit guiding — live</title>
<style>
  html, body { margin: 0; height: 100%; background: #111; color: #ccc;
    font: 13px -apple-system, Helvetica, Arial, sans-serif; }
  #bar { position: fixed; top: 0; left: 0; right: 0; z-index: 10;
    padding: 5px 10px; background: rgba(0,0,0,.65); }
  #bar b { color: #fff; }
  #img { display: block; max-width: 100vw; max-height: 100vh;
    margin: 0 auto; padding-top: 26px; box-sizing: border-box; }
</style>
</head>
<body>
<div id="bar">
  <b>slit guiding — live</b>
  &nbsp;·&nbsp; auto-refresh every <span id="iv">2</span>s
  &nbsp;·&nbsp; last refresh: <span id="t">—</span>
  &nbsp;·&nbsp; <a href="box_overlay.png" target="_blank" style="color:#6cf">box layout &#8599;</a>
</div>
<img id="img" src="dx_dy_vs_frame.png" alt="guide plot (waiting for first redraw)">
<script>
  var INTERVAL_MS = 2000;
  var SRC = "dx_dy_vs_frame.png";
  function refresh() {
    var next = new Image();
    next.onload = function () {
      document.getElementById("img").src = next.src;
      document.getElementById("t").textContent = new Date().toLocaleTimeString();
    };
    next.onerror = function () {};
    next.src = SRC + "?t=" + Date.now();
  }
  document.getElementById("iv").textContent = (INTERVAL_MS / 1000);
  setInterval(refresh, INTERVAL_MS);
  refresh();
</script>
</body>
</html>
"""


def open_in_browser(path):
    """Best-effort: pop the page up in the default browser (open / xdg-open).
    Never raises; returns True if an opener was launched."""
    import shutil
    import subprocess
    if sys.platform == "darwin":
        opener = "open"
    elif sys.platform.startswith("linux"):
        opener = shutil.which("xdg-open")
    else:
        opener = None
    if not opener:
        return False
    try:
        subprocess.Popen([opener, path],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception:
        return False


def ensure_live_html(out_dir):
    try:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, LIVE_HTML_NAME), "w") as fh:
            fh.write(LIVE_HTML_DOC)
    except OSError as e:
        print(f"  (could not write live.html: {e})")


def frame_fingerprint(cfg):
    paths = glob.glob(os.path.join(cfg["source"], "hen[0-9][0-9][0-9][0-9].fits"))
    frames = []
    for p in paths:
        try:
            fno = int(os.path.basename(p)[3:7])
        except ValueError:
            continue
        if fno >= cfg["start_frame"]:
            frames.append((fno, p))
    if not frames:
        return (0, -1, -1)
    max_fno, newest = max(frames, key=lambda t: t[0])
    try:
        size = os.path.getsize(newest)
    except OSError:
        size = -1
    return (len(frames), max_fno, size)


def _mtime(path):
    try:
        return os.path.getmtime(path)
    except OSError:
        return -1.0


def decide_reason(fp, ev, cf, last_frame, last_events, last_config):
    if last_frame is None:
        return "startup"
    reasons = []
    if fp != last_frame:
        reasons.append(f"new/updated frame (max {fp[1]})")
    if ev != last_events:
        reasons.append("motion_events.txt edit")
    if cf != last_config:
        reasons.append("config.txt edit")
    return " + ".join(reasons) if reasons else None


def _stamp():
    return time.strftime("%H:%M:%S")


def print_status(summary, cfg):
    if summary is None:
        return
    lo, hi = summary["frames"]
    print(f"  frames {lo}..{hi} (n={summary['n']})   "
          f"new:{summary['n_new']} cached:{summary['n'] - summary['n_new']}")
    t, c = summary["target"], summary["comp"]
    tdx, tdy = summary["latest"][t]
    cdx, cdy = summary["latest"][c]
    print(f"  {t:<12s} dx {tdx:+.2f}  dy {tdy:+.2f}        "
          f"{c:<12s} dx {cdx:+.2f}  dy {cdy:+.2f}")
    ddx, ddy = summary["diff_latest"]
    print(f"  diff{'':<8s} dx {ddx:+.2f}  dy {ddy:+.2f}   |d|={np.hypot(ddx, ddy):.2f} px")
    if np.isfinite(tdx) and np.isfinite(tdy):
        rec = keypress.recommend(tdx, tdy, cfg["pa"])
        print(f"  ==> PRESS:  {rec['text']}")


def watch(run_dir, *, last_n=30, interval=5.0, open_browser=True):
    # Line-buffer stdout so the status readout streams even when piped/redirected
    # (e.g. `watch | tee night.log`), not just when attached to a terminal.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    config_path = os.path.join(run_dir, "config.txt")
    events_file = os.path.join(run_dir, "motion_events.txt")
    out = os.path.join(run_dir, "outputs")
    ensure_live_html(out)
    live = os.path.join(out, LIVE_HTML_NAME)
    cfg0 = cfgmod.load_config(config_path)
    print(f"watch: {run_dir}")
    print(f"       source: {cfg0['source']}")
    print(f"       poll every {interval:g}s, --last-n {last_n}, PA={cfg0['pa']}")
    print(f"       VIEW (auto-updating):  open {live}")
    print(f"       (macOS Preview won't auto-refresh — use the browser page.)")
    if open_browser:
        opened = open_in_browser(live)
        print(f"       {'opened it in your browser.' if opened else 'auto-open unavailable; open it manually.'}")
    print(f"       Ctrl-C to stop.\n")
    if not os.path.isdir(cfg0["source"]):
        print(f"WARNING: source not found: {cfg0['source']} — set `source` in "
              f"config.txt (will keep polling).\n")

    last_frame = last_events = last_config = None
    try:
        while True:
            cfg = cfgmod.load_config(config_path)
            fp = frame_fingerprint(cfg)
            ev = _mtime(events_file)
            cf = _mtime(config_path)
            reason = decide_reason(fp, ev, cf, last_frame, last_events, last_config)
            if reason is not None:
                config_changed = (last_config is None) or (cf != last_config)
                last_frame, last_events, last_config = fp, ev, cf
                print(f"\n{'=' * 64}\n{_stamp()}  redraw — {reason}")
                try:
                    if config_changed:
                        montage.write_box_overlay(run_dir, cfg)
                    summary = montage.run_montage(run_dir, cfg, last_n=last_n)
                    print_status(summary, cfg)
                except Exception as e:  # noqa: BLE001 - keep the loop alive
                    print(f"  redraw failed: {e}")
            time.sleep(interval)
    except KeyboardInterrupt:
        print(f"\n{_stamp()}  watch stopped.")
