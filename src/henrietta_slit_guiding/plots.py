"""
Post-night plots driven by motion_events.txt:
  path  — quiver trail of the cumulative keypress nudges
  rate  — histogram of commanded motion (arcsec) vs time
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.time import Time

from . import config as cfgmod

ALIASES = {
    "u": "up", "n": "up", "north": "up", "d": "down", "s": "down", "south": "down",
    "r": "right", "e": "right", "east": "right", "l": "left", "w": "left", "west": "left",
}


def _canon(d):
    return ALIASES.get(d.lower(), d.lower())


def _load_events(events_file):
    events = []
    if not os.path.exists(events_file):
        return events
    with open(events_file) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 3)
            if len(parts) < 3:
                continue
            try:
                frame, amount = int(parts[0]), float(parts[2])
            except ValueError:
                continue
            events.append(dict(frame=frame, direction=_canon(parts[1]), amount=amount))
    return events


def _step(direction, amount):
    if direction == "up":
        return 0.0, -amount
    if direction == "down":
        return 0.0, +amount
    if direction == "right":
        return +amount, 0.0
    if direction == "left":
        return -amount, 0.0
    raise ValueError(f"unknown direction {direction!r}")


def _mjd_for_frame(cfg, fno):
    p = os.path.join(cfg["source"], f"hen{fno:04d}.fits")
    if not os.path.exists(p):
        return float("nan")
    try:
        return Time(str(fits.getheader(p).get("DATE-OBS", ""))).mjd
    except Exception:
        return float("nan")


def _tag(cfg):
    t, c = cfgmod.target_comp(cfg["boxes"])
    return f"{t} + {c}"


def plot_path(run_dir, cfg):
    out = os.path.join(run_dir, "outputs")
    os.makedirs(out, exist_ok=True)
    events = _load_events(os.path.join(run_dir, "motion_events.txt"))
    if not events:
        print("no motion events to plot"); return
    pixscale, pa = cfg["pixscale"], cfg["pa"]
    xs, ys, dirs = [0.0], [0.0], [None]
    for ev in events:
        ddx, ddy = _step(ev["direction"], ev["amount"])
        xs.append(xs[-1] + ddx); ys.append(ys[-1] + ddy); dirs.append(ev["direction"])
    xs, ys = np.array(xs), np.array(ys)
    u, v = np.diff(xs), np.diff(ys)
    cmap = {"up": "C0", "down": "C1", "left": "C2", "right": "C3"}
    cols = [cmap.get(d, "k") for d in dirs[1:]]
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.quiver(xs[:-1], ys[:-1], u, v, color=cols, angles="xy",
              scale_units="xy", scale=1.0, width=0.005, headwidth=4.5)
    ax.plot(xs[0], ys[0], "o", color="k", ms=10, label=f"start (0,0)")
    ax.plot(xs[-1], ys[-1], "*", color="gold", ms=18, mec="k",
            label=f"end ({xs[-1]:+.0f},{ys[-1]:+.0f})")
    ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_xlabel("guider X [pix, rel.]"); ax.set_ylabel("guider Y [pix, rel.]")
    net = (xs[-1] - xs[0], ys[-1] - ys[0])
    ax.set_title(f"{_tag(cfg)} — keypress trail   "
                 f"PA={pa:.2f}°   {pixscale}\"/pix\n"
                 f"N={len(events)} events   "
                 f"net Δ=({net[0]:+.0f},{net[1]:+.0f}) px = "
                 f"({net[0]*pixscale:+.2f},{net[1]*pixscale:+.2f})\"")
    fig.tight_layout()
    p = os.path.join(out, "guider_path.png")
    fig.savefig(p, dpi=120); plt.close(fig)
    print("wrote", p)


def plot_rate(run_dir, cfg, bin_minutes=20):
    out = os.path.join(run_dir, "outputs")
    os.makedirs(out, exist_ok=True)
    events = _load_events(os.path.join(run_dir, "motion_events.txt"))
    if not events:
        print("no motion events to plot"); return
    pixscale = cfg["pixscale"]
    mjds = np.array([_mjd_for_frame(cfg, ev["frame"]) for ev in events])
    dirs = np.array([ev["direction"] for ev in events])
    amts = np.array([ev["amount"] for ev in events])
    ok = np.isfinite(mjds)
    if not ok.any():
        print("no usable frame timestamps"); return
    mjds, dirs, amts = mjds[ok], dirs[ok], amts[ok]
    t_min = (mjds - mjds.min()) * 1440.0
    t0 = Time(mjds.min(), format="mjd")
    n_bins = max(1, int(np.ceil(float(t_min.max()) / bin_minutes)))
    edges = np.arange(n_bins + 1) * bin_minutes
    centers = 0.5 * (edges[:-1] + edges[1:])
    is_x = (dirs == "left") | (dirs == "right")
    is_y = (dirs == "up") | (dirs == "down")
    motion = np.abs(amts) * pixscale
    hx, _ = np.histogram(t_min[is_x], bins=edges, weights=motion[is_x])
    hy, _ = np.histogram(t_min[is_y], bins=edges, weights=motion[is_y])
    htot, _ = np.histogram(t_min, bins=edges, weights=motion)
    fig, ax = plt.subplots(figsize=(11, 5))
    w = bin_minutes * 0.27
    ax.bar(centers - w, hx, width=w, label="X (L+R)", color="C2", edgecolor="k", lw=0.4)
    ax.bar(centers, hy, width=w, label="Y (U+D)", color="C0", edgecolor="k", lw=0.4)
    ax.bar(centers + w, htot, width=w, label="total", color="0.4", edgecolor="k", lw=0.4, alpha=0.7)
    ax.set_xticks(centers)
    ax.set_xticklabels([(t0 + c / 1440.0).iso[11:16] for c in centers], rotation=45, ha="right", fontsize=8)
    ax.set_xlabel(f"UT (bins {bin_minutes} min; first at {t0.iso[11:19]})")
    ax.set_ylabel(f"commanded motion [arcsec] / {bin_minutes} min")
    ax.grid(axis="y", alpha=0.3); ax.legend()
    ax.set_title(f"{_tag(cfg)} — commanded motion vs time   "
                 f"sum X={hx.sum():.2f}\" Y={hy.sum():.2f}\" total={htot.sum():.2f}\"")
    fig.tight_layout()
    p = os.path.join(out, "motion_rate.png")
    fig.savefig(p, dpi=120); plt.close(fig)
    print("wrote", p)
