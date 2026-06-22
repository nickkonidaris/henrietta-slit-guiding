"""
Dump headers of the frames in range, to pick `start_frame` and confirm the
`object` / `filter` strings before editing config.txt.
"""
from __future__ import annotations

import glob
import os

from astropy.io import fits


def dump_headers(cfg: dict) -> None:
    src = cfg["source"]
    f0 = cfg["start_frame"]
    paths = sorted(glob.glob(os.path.join(src, "hen[0-9][0-9][0-9][0-9].fits")))
    rows = []
    for p in paths:
        try:
            f = int(os.path.basename(p)[3:7])
        except ValueError:
            continue
        if f >= f0:
            rows.append((f, p))
    if not rows:
        print(f"no hen####.fits at/after frame {f0} in {src!r} "
              f"(set `source`/`start_frame` in config.txt).")
        return
    print(f"{len(rows)} frames at/after {f0} in {src}")
    print(f"{'frame':>5s} {'OBJECT':<40s} {'FILTER':<6s} {'EXP':>6s} "
          f"{'AM':>5s} {'ROT':>8s} {'DATE-OBS':<23s}")
    for f, p in rows:
        h = fits.getheader(p)
        print(f"{f:>5d} {str(h.get('OBJECT', '?'))[:40]:<40s} "
              f"{str(h.get('FILTER', '?'))[:6]:<6s} "
              f"{h.get('EXPTIME', 0):>6.1f} {h.get('AIRMASS', 0):>5.3f} "
              f"{h.get('ROTANGLE', 0):>8.3f} {str(h.get('DATE-OBS', '?'))[:23]:<23s}")
