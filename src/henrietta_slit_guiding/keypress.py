"""
Translate a measured detector drift (dx, dy) into TCS arrow keypresses that
re-center the star, at the night's instrument PA.

The WASP-87 calibration gives each key's per-unit detector response at
PA0=79.46 deg.  TCS commands are sky-fixed, so the response rotates with PA
by theta = PA0 - PA.  Solving M·n = -drift gives the integer up/left presses
(negative -> down/right).
"""
from __future__ import annotations

import os

import numpy as np

# --- WASP-87 calibration (merge-opposite), px per TCS unit at PA0 ---
PA0 = 79.46
V0_UP = (0.394, -0.280)
V0_LEFT = (0.122, 0.509)


def rotate(v, theta_rad):
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    vx, vy = v
    return (vx * c - vy * s, vx * s + vy * c)


def response_matrix(pa, scale=1.0):
    theta = np.deg2rad(PA0 - pa)
    v_up = rotate(V0_UP, theta)
    v_left = rotate(V0_LEFT, theta)
    M = np.array([[v_up[0], v_left[0]], [v_up[1], v_left[1]]]) * scale
    return M, v_up, v_left


def recommend(dx, dy, pa, scale=1.0, deadband=0.2):
    """Return a dict describing the correction for drift (dx, dy) at this PA.

    Keys: presses (list like ['UP 1','LEFT 1']), text (one-line summary),
    within_deadband, dmag, and the exact n_up/n_left.
    """
    dmag = float(np.hypot(dx, dy))
    if dmag < deadband:
        return dict(within_deadband=True, presses=[], dmag=dmag,
                    text=f"within deadband ({deadband}px) — hold")
    M, _, _ = response_matrix(pa, scale)
    n = np.linalg.solve(M, np.array([-dx, -dy]))
    n_up, n_left = float(n[0]), float(n[1])
    ud_dir = "up" if n_up >= 0 else "down"
    lr_dir = "left" if n_left >= 0 else "right"
    ud_r, lr_r = int(round(abs(n_up))), int(round(abs(n_left)))
    presses = []
    if ud_r > 0:
        presses.append(f"{ud_dir.upper()} {ud_r}")
    if lr_r > 0:
        presses.append(f"{lr_dir.upper()} {lr_r}")
    text = "   ".join(presses) if presses else "sub-unit nudge — hold"
    return dict(within_deadband=False, presses=presses, dmag=dmag, text=text,
                n_up=n_up, n_left=n_left, ud_dir=ud_dir, lr_dir=lr_dir)


def latest_drift(npz_path, star):
    """Most recent finite (dx, dy, frame) for `star` from the npz, or None."""
    if not os.path.exists(npz_path):
        return None
    z = np.load(npz_path)
    kx, ky = f"dx_{star}", f"dy_{star}"
    if kx not in z or ky not in z:
        raise SystemExit(f"{kx}/{ky} not in {npz_path}; keys: {list(z.keys())}")
    dx, dy, frames = z[kx], z[ky], z["frames"]
    good = np.where(np.isfinite(dx) & np.isfinite(dy))[0]
    if len(good) == 0:
        return None
    i = good[-1]
    return float(dx[i]), float(dy[i]), int(frames[i])


def run_keypress(run_dir, cfg, *, dx=None, dy=None, scale=1.0, deadband=0.2,
                 star=None, pa=None):
    """CLI entry: print the recommended presses for an explicit or latest drift."""
    from . import config as cfgmod
    from .montage import NPZ_NAME

    pa = cfg["pa"] if pa is None else pa
    if star is None:
        star, _ = cfgmod.target_comp(cfg["boxes"])

    if dx is not None and dy is not None:
        src = "given"
        frame = None
    else:
        npz = os.path.join(run_dir, "outputs", NPZ_NAME)
        got = latest_drift(npz, star)
        if got is None:
            raise SystemExit(f"no --dx/--dy and no finite drift in {npz}; "
                             f"run `watch` first or pass --dx/--dy.")
        dx, dy, frame = got
        src = f"npz latest ({star}, frame {frame})"

    M, v_up, v_left = response_matrix(pa, scale)
    print(f"instrument PA = {pa:.2f} deg   (theta = PA0 - PA = {PA0 - pa:.2f})   scale = {scale}")
    print(f"per-unit detector response [px/unit]:")
    print(f"   UP   -> (dx {v_up[0]*scale:+.3f}, dy {v_up[1]*scale:+.3f})")
    print(f"   LEFT -> (dx {v_left[0]*scale:+.3f}, dy {v_left[1]*scale:+.3f})")
    print(f"\nmeasured drift [{src}]:  dx = {dx:+.3f}  dy = {dy:+.3f}  "
          f"(|d| = {np.hypot(dx, dy):.3f} px)")
    rec = recommend(dx, dy, pa, scale, deadband)
    if rec["within_deadband"]:
        print(f"\n--> {rec['text']}: NO correction needed.")
        return
    print(f"\nexact correction:  {rec['ud_dir']} {abs(rec['n_up']):.2f}   "
          f"{rec['lr_dir']} {abs(rec['n_left']):.2f}  (TCS units)")
    if rec["presses"]:
        print(f"\n==> PRESS:  {'   '.join(rec['presses'])}")
    else:
        print("\n--> rounds to 0 presses; sub-unit nudge, hold.")
