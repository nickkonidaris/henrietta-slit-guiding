"""
Find the two trace columns in the first matching frame, print suggested
`box` lines (in config.txt format), and save a figure showing HOW they were
found: the Y-smashed column profile with the detected peaks and the chosen
extraction boxes marked.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std, sigma_clip
from scipy.signal import find_peaks

from . import config as cfgmod
from .montage import discover_frames

HW = 21


def find_stars(cfg: dict, run_dir: str | None = None) -> None:
    boxes = cfg["boxes"]
    target, comp = cfgmod.target_comp(boxes)
    y_lo = boxes[target]["y_lo"]
    y_hi = boxes[target]["y_hi"]

    avail = discover_frames(cfg)
    if not avail:
        print(f"no matching frames (object={cfg['object']!r} filter={cfg['filter']!r} "
              f"start_frame={cfg['start_frame']} source={cfg['source']!r}) — "
              f"set those in config.txt first.")
        return
    fno = avail[0]
    path = os.path.join(cfg["source"], f"hen{fno:04d}.fits")
    bpm = fits.getdata(cfg["bpm_path"]).astype(bool)
    img = np.where(bpm, fits.getdata(path).astype("f4"), np.nan)
    hdr = fits.getheader(path)
    print(f"frame {fno}: OBJECT={hdr.get('OBJECT')!r} FILTER={hdr.get('FILTER')!r} "
          f"ROTANGLE={hdr.get('ROTANGLE')}")
    print(f"  Y range = {y_lo}..{y_hi}")

    sub = img[y_lo:y_hi, :]
    masked = sigma_clip(sub, sigma=5.0, axis=0, masked=True,
                        cenfunc="median", stdfunc=mad_std, maxiters=2)
    x_prof = (np.ma.mean(masked, axis=0).filled(0.0) * sub.shape[0]).astype("f8")
    x_prof_pos = np.where(x_prof > 0, x_prof, 0.0)
    peaks, _ = find_peaks(x_prof_pos, distance=80,
                          prominence=np.nanpercentile(x_prof_pos, 95) * 0.3)
    if len(peaks) < 2:
        print(f"found < 2 peaks ({len(peaks)}); adjust the Y range in config.txt.")
        return
    order = np.argsort(x_prof_pos[peaks])[::-1]
    top = peaks[order][:5]
    bright_col, faint_col = int(top[0]), int(top[1])
    print("\nTop peak columns (brightest first):")
    for c in top:
        print(f"  col {int(c):>5d}   value {x_prof_pos[c]:>10.0f}")
    print(f"\n  brighter trace: col {bright_col}   fainter: col {faint_col}")
    print("\nPaste into config.txt (brightness guess; swap if the field is mirrored):")
    print("  #    name        role    x_center  x_halfwidth  y_lo  y_hi")
    print(f"  box  {target:<10s} target  {faint_col:<8d}  {HW:<11d}  {y_lo:<4d}  {y_hi}")
    print(f"  box  {comp:<10s} comp    {bright_col:<8d}  {HW:<11d}  {y_lo:<4d}  {y_hi}")
    print("  # (target assumed fainter; verify against the field)")

    # --- figure: how the peaks were found ---
    if run_dir is None:
        run_dir = os.getcwd()
    out = os.path.join(run_dir, "outputs")
    os.makedirs(out, exist_ok=True)
    cols = np.arange(len(x_prof_pos))
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(cols, x_prof_pos, "-", color="0.4", lw=0.8,
            label=f"Y-smashed profile (rows {y_lo}..{y_hi})")
    ax.plot(peaks, x_prof_pos[peaks], "x", color="C7", ms=7, label="detected peaks")
    assign = [(faint_col, target, "target", "C3"), (bright_col, comp, "comp", "C0")]
    for col, name, role, color in assign:
        ax.axvspan(col - HW, col + HW + 1, color=color, alpha=0.18)
        ax.axvline(col, color=color, lw=1.2)
        ax.annotate(f"{name}\n({role}) col {col}",
                    xy=(col, x_prof_pos[col]), xytext=(0, 12),
                    textcoords="offset points", ha="center", va="bottom",
                    color=color, fontsize=9, weight="bold")
    ax.set_xlabel("detector column [pix]")
    ax.set_ylabel("Y-smashed signal")
    ax.set_title(f"find-stars: frame {fno}  OBJECT={hdr.get('OBJECT', '?')}  "
                 f"(shaded = suggested boxes, ±{HW} px)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper center", fontsize=8)
    fig.tight_layout()
    p = os.path.join(out, "find_stars.png")
    fig.savefig(p, dpi=120)
    plt.close(fig)
    print(f"\nwrote {p}")
