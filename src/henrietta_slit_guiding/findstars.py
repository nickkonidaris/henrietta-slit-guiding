"""
Find the two trace columns in the first matching frame and print suggested
`box` lines (in config.txt format) to paste in.  Uses the target/comp names
already in the config so the output is directly pasteable.
"""
from __future__ import annotations

import numpy as np
from astropy.io import fits
from astropy.stats import mad_std, sigma_clip
from scipy.signal import find_peaks

from . import config as cfgmod
from .montage import discover_frames, HP_PRE_FIT  # noqa: F401  (HP unused; keeps import parity)

HW = 21


def find_stars(cfg: dict) -> None:
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
    import os
    path = os.path.join(cfg["source"], f"hen{fno:04d}.fits")
    bpm = fits.getdata(cfg["bpm_path"]).astype(bool)
    img = fits.getdata(path).astype("f4")
    img = np.where(bpm, img, np.nan)
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
    print("\nTop peak columns (brightest first):")
    for c in top:
        print(f"  col {int(c):>5d}   value {x_prof_pos[c]:>10.0f}")

    bright_col, faint_col = int(top[0]), int(top[1])
    print(f"\n  brighter trace: col {bright_col}   fainter: col {faint_col}")
    print("\nPaste into config.txt (brightness guess; swap if the field is mirrored):")
    print("  #    name        role    x_center  x_halfwidth  y_lo  y_hi")
    print(f"  box  {target:<10s} target  {faint_col:<8d}  {HW:<11d}  {y_lo:<4d}  {y_hi}")
    print(f"  box  {comp:<10s} comp    {bright_col:<8d}  {HW:<11d}  {y_lo:<4d}  {y_hi}")
    print("  # (target assumed fainter; verify against the field)")
