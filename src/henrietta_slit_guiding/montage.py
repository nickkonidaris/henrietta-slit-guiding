"""
Per-frame slit-motion analysis: X drift via 1-D cross-correlation of the
Y-smashed column profile against a template, Y drift via an inverted-Gaussian
fit to an absorption dip in the spectrum.  Differential (target - comp) is
the photometrically relevant signal.

Only NEW frames are fit; results are cached per frame in
outputs/motion_cache.json (reused when the file size + a signature match).
"""
from __future__ import annotations

import json
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std, sigma_clip
from astropy.time import Time
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter, uniform_filter1d

from . import config as cfgmod

# --- analysis constants (instrument/method, not per-night) ---
F1 = 9999
SEARCH = 8
SIGMA_CLIP_N = 5.0
TOPHAT = 3
GAUSS_GUESS_ROW = 40
GAUSS_HALF_WINDOW = 15
SIGMA_MAX = 8.0
HP_PRE_FIT = 41
N_COLS = 12
EVENT_FRAME_OFFSET = -1
REQUIRE_OBSERVATIONS = False
ARROWS = {
    "up": "↑", "down": "↓", "left": "←", "right": "→",
    "u": "↑", "d": "↓", "l": "←", "r": "→",
    "n": "↑", "s": "↓", "e": "→", "w": "←",
    "north": "↑", "south": "↓", "east": "→", "west": "←",
}

PNG_NAME = "dx_dy_vs_frame.png"
NPZ_NAME = "motion_montage.npz"
CACHE_NAME = "motion_cache.json"
OVERLAY_NAME = "box_overlay.png"


# ---------- frame access ----------
def slope_path(cfg: dict, fno: int) -> str:
    return os.path.join(cfg["source"], f"hen{fno:04d}.fits")


def load_signal(cfg: dict, fno: int):
    p = slope_path(cfg, fno)
    return fits.getdata(p).astype("f4"), fits.getheader(p)


def is_target_obs(cfg: dict, hdr: fits.Header) -> bool:
    obj = str(hdr.get("OBJECT", ""))
    if cfg["object"] and cfg["object"].lower() not in obj.lower():
        return False
    if REQUIRE_OBSERVATIONS and "Observations" not in obj:
        return False
    return True


def discover_frames(cfg: dict) -> list:
    avail = []
    want = cfg["filter"]
    for f in range(cfg["start_frame"], F1 + 1):
        p = slope_path(cfg, f)
        if not os.path.exists(p):
            continue
        try:
            h = fits.getheader(p)
        except Exception:
            continue
        if is_target_obs(cfg, h) and (not want or h.get("FILTER") == want):
            avail.append(f)
    return avail


# ---------- profiles / fits (pure) ----------
def extract_box(img, good_full, box):
    s = img[box["y_lo"]:box["y_hi"], box["x_lo"]:box["x_hi"]].copy()
    g = good_full[box["y_lo"]:box["y_hi"], box["x_lo"]:box["x_hi"]].copy()
    s = np.where(np.isfinite(s) & g, s, 0.0)
    n_x = s.shape[1]
    edge = max(1, n_x // 6)
    sky = np.nanmedian(np.concatenate([s[:, :edge], s[:, -edge:]], axis=1),
                       axis=1, keepdims=True)
    sky = np.where(np.isfinite(sky), sky, 0.0)
    return s - sky, g, s, sky


def make_x_profile(sub):
    n_y = sub.shape[0]
    masked = sigma_clip(sub, sigma=SIGMA_CLIP_N, axis=0, masked=True,
                        cenfunc="median", stdfunc=mad_std, maxiters=2)
    per_col = np.ma.mean(masked, axis=0).filled(0.0)
    return (per_col * n_y).astype("f8")


def make_y_spec(sub):
    n_x = sub.shape[1]
    masked = sigma_clip(sub, sigma=SIGMA_CLIP_N, axis=1, masked=True,
                        cenfunc="median", stdfunc="std", maxiters=2)
    per_row = np.ma.mean(masked, axis=1).filled(0.0)
    return (per_row * n_x).astype("f8")


def tophat(spec, n: int = TOPHAT):
    if n <= 1:
        return spec.copy()
    return uniform_filter1d(spec, size=n, mode="nearest")


def gaussian_dip(x, c, A, mu, sigma):
    return c + A * np.exp(-((x - mu) ** 2) / (2.0 * sigma ** 2))


def hp_continuum_subtract(spec, npx: int = HP_PRE_FIT):
    if npx <= 1:
        return spec.copy()
    cont = median_filter(np.nan_to_num(spec, nan=np.nanmedian(spec)),
                         size=npx, mode="reflect")
    return spec - cont


def fit_dip(spec, guess_row: int = GAUSS_GUESS_ROW, half_window: int = GAUSS_HALF_WINDOW):
    n = len(spec)
    lo = max(0, guess_row - half_window)
    hi = min(n, guess_row + half_window + 1)
    flat = hp_continuum_subtract(spec)
    x = np.arange(lo, hi, dtype="f8")
    y = flat[lo:hi].astype("f8")
    if not np.all(np.isfinite(y)):
        m = np.isfinite(y)
        x = x[m]; y = y[m]
        if len(x) < 6:
            raise RuntimeError("not enough finite samples")
    c0 = 0.0
    i_min = int(np.argmin(y))
    A0 = float(y[i_min] - c0)
    mu0 = float(x[i_min])
    p0 = [c0, A0, mu0, 2.5]
    try:
        popt, _ = curve_fit(gaussian_dip, x, y, p0=p0,
                            bounds=([-np.inf, -np.inf, lo, 0.5],
                                    [np.inf, 0.0, hi, SIGMA_MAX]),
                            maxfev=4000)
    except Exception as e:
        raise RuntimeError(f"curve_fit failed: {e}")
    c, A, mu, sigma = (float(p) for p in popt)
    if not (lo <= mu <= hi):
        raise RuntimeError(f"mu={mu:.2f} out of fit window [{lo},{hi}]")
    return c, A, mu, sigma, x


def xcor_1d(profile, template, search: int = SEARCH):
    n = len(profile)
    lags = np.arange(-search, search + 1)
    C = np.zeros(len(lags))
    for k, L in enumerate(lags):
        if L >= 0:
            C[k] = float(np.sum(template[:n - L] * profile[L:]))
        else:
            C[k] = float(np.sum(template[-L:] * profile[:n + L]))
    k = int(np.argmax(C))
    if k == 0 or k == len(C) - 1:
        return float(lags[k])
    a, b, c = C[k - 1], C[k], C[k + 1]
    denom = a - 2.0 * b + c
    sub = 0.5 * (a - c) / denom if denom != 0 else 0.0
    return float(lags[k] + sub)


# ---------- events ----------
def load_events(events_file: str) -> list:
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
                frame = int(parts[0])
            except ValueError:
                continue
            events.append(dict(frame=frame, direction=parts[1], amount=parts[2],
                               comment=parts[3] if len(parts) >= 4 else ""))
    return events


def mark_events(ax, events):
    if not events:
        return
    ymin, ymax = ax.get_ylim()
    label_y = ymax - 0.06 * (ymax - ymin)
    by_x = {}
    for ev in events:
        x = ev["frame"] + EVENT_FRAME_OFFSET
        arrow = ARROWS.get(ev["direction"].lower(), ev["direction"])
        by_x.setdefault(x, []).append(f"{arrow}{ev['amount']}")
    for x, labels in by_x.items():
        ax.axvline(x, color="magenta", lw=0.8, ls="--", alpha=0.6)
        ax.text(x, label_y, "\n".join(labels), ha="center", va="top",
                color="magenta", fontsize=8, clip_on=True,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          alpha=0.8, edgecolor="magenta", lw=0.5))


# ---------- one-frame fit ----------
def fit_frame(cfg, fno, bpm, templates, template_mu, boxes):
    try:
        img, hdr = load_signal(cfg, fno)
    except Exception as e:
        print(f"  frame {fno} skip: {e}")
        return None
    img = np.where(bpm, img, np.nan)
    try:
        t_mjd = float(Time(hdr["DATE-OBS"]).mjd)
    except Exception:
        t_mjd = float("nan")
    stars, extras = {}, {}
    for star, box in boxes.items():
        sub_d, _, raw, sky_pr = extract_box(img, bpm, box)
        x_p = make_x_profile(sub_d)
        y_smooth = tophat(make_y_spec(sub_d))
        dxv = float(xcor_1d(x_p - x_p.mean(), templates[star]["x_t"]))
        muv = sigv = Av = cv = float("nan")
        dyv = float("nan")
        try:
            cv, Av, muv, sigv, _ = fit_dip(y_smooth)
            dyv = muv - template_mu[star]
        except RuntimeError as e:
            print(f"  frame {fno} {star} fit fail: {e}")
        stars[star] = dict(dx=dxv, dy=float(dyv), mu=float(muv), sigma=float(sigv),
                           A=float(Av), c=float(cv), obj=float(np.nansum(sub_d)),
                           sky=float(np.nanmedian(sky_pr)))
        extras[star] = dict(raw=raw, y_smooth=y_smooth)
    try:
        size = int(os.path.getsize(slope_path(cfg, fno)))
    except OSError:
        size = -1
    return dict(time_mjd=t_mjd, size=size, stars=stars, extras=extras)


# ---------- cache ----------
def cache_signature(cfg, f_template, boxes):
    return {
        "template_frame": int(f_template),
        "boxes": {s: [boxes[s]["x_lo"], boxes[s]["x_hi"], boxes[s]["y_lo"],
                      boxes[s]["y_hi"]] for s in boxes},
        "params": [SEARCH, float(SIGMA_CLIP_N), TOPHAT, GAUSS_GUESS_ROW,
                   GAUSS_HALF_WINDOW, HP_PRE_FIT, float(SIGMA_MAX)],
        "bpm": cfg["bpm_path"],
        "ver": 1,
    }


def load_cache(signature, path):
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as fh:
            blob = json.load(fh)
        if blob.get("signature") != signature:
            return {}
        return blob.get("frames", {})
    except Exception as e:
        print(f"  (cache ignored: {e})")
        return {}


def save_cache(signature, frames, path):
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump({"signature": signature, "frames": frames}, fh)
        os.replace(tmp, path)
    except Exception as e:
        print(f"  (could not write cache: {e})")


# ---------- box overlay ----------
def write_box_overlay(run_dir, cfg, out_path=None):
    boxes = cfg["boxes"]
    target, comp = cfgmod.target_comp(boxes)
    colors = {target: "C3", comp: "C0"}
    avail = discover_frames(cfg)
    if not avail:
        print(f"box overlay: no matching frames found (start_frame={cfg['start_frame']})")
        return None
    fno = avail[0]
    img, hdr = load_signal(cfg, fno)
    bpm = fits.getdata(cfg["bpm_path"]).astype(bool)
    img = np.where(bpm, img, np.nan)
    ny, nx = img.shape
    vmin = float(np.nanpercentile(img, 30))
    vmax = float(np.nanpercentile(img, 99.5))
    fig_w = 14.0
    fig_h = max(4.0, min(11.0, fig_w * ny / nx))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(img, origin="lower", cmap="gray", vmin=vmin, vmax=vmax,
              interpolation="nearest", aspect="auto")
    for star, box in boxes.items():
        col = colors.get(star, "C1")
        w = box["x_hi"] - box["x_lo"]
        h = box["y_hi"] - box["y_lo"]
        ax.add_patch(mpatches.Rectangle((box["x_lo"] - 0.5, box["y_lo"] - 0.5),
                                        w, h, fill=False, edgecolor=col, lw=1.8))
        ax.text(0.5 * (box["x_lo"] + box["x_hi"]), min(ny - 2, box["y_hi"] + 8),
                f"{star} ({box.get('role', '')})\n"
                f"x {box['x_lo']}..{box['x_hi']}  y {box['y_lo']}..{box['y_hi']}",
                color=col, ha="center", va="bottom", fontsize=9, weight="bold")
    ax.set_xlim(0, nx); ax.set_ylim(0, ny)
    ax.set_xlabel("detector X [pix]"); ax.set_ylabel("detector Y [pix]")
    ax.set_title(f"Extraction boxes on frame {fno}   "
                 f"OBJECT={hdr.get('OBJECT', '?')}  FILTER={hdr.get('FILTER', '?')}\n"
                 f"(edit config.txt to move a box)")
    fig.tight_layout()
    out = os.path.join(run_dir, "outputs")
    os.makedirs(out, exist_ok=True)
    if out_path is None:
        out_path = os.path.join(out, OVERLAY_NAME)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"wrote {out_path}")
    return out_path


def _latest_finite(dx, dy):
    good = np.where(np.isfinite(dx) & np.isfinite(dy))[0]
    if len(good) == 0:
        return float("nan"), float("nan")
    i = good[-1]
    return float(dx[i]), float(dy[i])


# ---------- the main analysis ----------
def run_montage(run_dir, cfg, *, last_n=30, show_all=False, all_plots=False,
                no_cache=False):
    out = os.path.join(run_dir, "outputs")
    os.makedirs(out, exist_ok=True)
    events_file = os.path.join(run_dir, "motion_events.txt")
    boxes = cfg["boxes"]
    target, comp = cfgmod.target_comp(boxes)
    colors = {target: "C3", comp: "C0"}
    diff_name = f"{target} - {comp}"

    bpm = fits.getdata(cfg["bpm_path"]).astype(bool)
    avail = discover_frames(cfg)
    if not avail:
        print(f"no matching frames (object={cfg['object']!r} filter={cfg['filter']!r} "
              f"start_frame={cfg['start_frame']} source={cfg['source']!r})")
        return None
    print(f"{len(avail)} frames found: {avail[0]}..{avail[-1]}")

    f_template = avail[0]
    img0, _ = load_signal(cfg, f_template)
    img0 = np.where(bpm, img0, np.nan)
    templates, template_mu = {}, {}
    for star, box in boxes.items():
        sub_t, _, _, _ = extract_box(img0, bpm, box)
        x_t = make_x_profile(sub_t)
        y_t_smooth = tophat(make_y_spec(sub_t))
        try:
            c, A, mu, sigma, fit_x = fit_dip(y_t_smooth)
        except RuntimeError as e:
            print(f"  template {star}: dip fit FAILED: {e}")
            raise
        templates[star] = dict(x_t=x_t - x_t.mean(),
                               y_t_flat=hp_continuum_subtract(y_t_smooth),
                               mu=mu, sigma=sigma, A=A, c=c, fit_x=fit_x)
        template_mu[star] = mu
        print(f"  template {star:<12s} mu={mu:.3f} sigma={sigma:.2f}")

    nf = len(avail)
    dx = {s: np.full(nf, np.nan) for s in boxes}
    dy = {s: np.full(nf, np.nan) for s in boxes}
    mu_all = {s: np.full(nf, np.nan) for s in boxes}
    sigma_all = {s: np.full(nf, np.nan) for s in boxes}
    A_all = {s: np.full(nf, np.nan) for s in boxes}
    c_all = {s: np.full(nf, np.nan) for s in boxes}
    obj_counts = {s: np.full(nf, np.nan) for s in boxes}
    sky_counts = {s: np.full(nf, np.nan) for s in boxes}
    box_imgs = {s: [None] * nf for s in boxes}
    y_specs_smooth = {s: [None] * nf for s in boxes}
    times_mjd = np.full(nf, np.nan)

    use_cache = not all_plots and not no_cache
    sig = cache_signature(cfg, f_template, boxes)
    cache_file = os.path.join(out, CACHE_NAME)
    cache = load_cache(sig, cache_file) if use_cache else {}
    n_new = 0
    t0 = time.time()
    for i, fno in enumerate(avail):
        key = str(fno)
        entry = None
        if use_cache:
            cached = cache.get(key)
            if cached is not None:
                try:
                    cur_size = int(os.path.getsize(slope_path(cfg, fno)))
                except OSError:
                    cur_size = -1
                if cached.get("size") == cur_size:
                    entry = cached
        extras = None
        if entry is None:
            res = fit_frame(cfg, fno, bpm, templates, template_mu, boxes)
            if res is None:
                continue
            n_new += 1
            extras = res.get("extras")
            entry = res
            if use_cache:
                cache[key] = dict(time_mjd=res["time_mjd"], size=res["size"],
                                  stars=res["stars"])
        times_mjd[i] = entry["time_mjd"]
        for star in boxes:
            s = entry["stars"][star]
            dx[star][i] = s["dx"]; dy[star][i] = s["dy"]
            mu_all[star][i] = s["mu"]; sigma_all[star][i] = s["sigma"]
            A_all[star][i] = s["A"]; c_all[star][i] = s["c"]
            obj_counts[star][i] = s["obj"]; sky_counts[star][i] = s["sky"]
            if extras is not None:
                box_imgs[star][i] = extras[star]["raw"]
                y_specs_smooth[star][i] = extras[star]["y_smooth"]
        if i % 20 == 0 or i == nf - 1:
            print(f"  frame {fno} ({i+1}/{nf})  t={time.time()-t0:.1f}s")
    if use_cache:
        save_cache(sig, cache, cache_file)
    print(f"  ({n_new} new frames fit, {nf - n_new} reused from cache)")

    dxd = dx[target] - dx[comp]
    dyd = dy[target] - dy[comp]
    print("\n--- per-star summary ---")
    for star in boxes:
        print(f"  {star:<12s} dx σ={np.nanstd(dx[star]):.4f}  "
              f"dy σ={np.nanstd(dy[star]):.4f}")
    print(f"  diff ({diff_name})  dx σ={np.nanstd(dxd):.4f}  dy σ={np.nanstd(dyd):.4f}")

    _plot_timeseries(out, avail, dx, dy, dxd, dyd, obj_counts, sky_counts,
                     colors, diff_name, boxes, f_template, events_file,
                     last_n, show_all)
    if all_plots:
        _plot_box_montages(out, avail, boxes, box_imgs, dx, dy, mu_all,
                           template_mu, f_template)

    per_star = {}
    for s in boxes:
        per_star[f"dx_{s}"] = dx[s]; per_star[f"dy_{s}"] = dy[s]
        per_star[f"mu_{s}"] = mu_all[s]; per_star[f"sigma_{s}"] = sigma_all[s]
        per_star[f"obj_{s}"] = obj_counts[s]; per_star[f"sky_{s}"] = sky_counts[s]
    npz_path = os.path.join(out, NPZ_NAME)
    np.savez(npz_path, frames=np.array(avail), template_frame=f_template,
             times_mjd=times_mjd, names=np.array(list(boxes)),
             target=target, comp=comp,
             template_mu=np.array([template_mu[s] for s in boxes]), **per_star)

    latest = {s: _latest_finite(dx[s], dy[s]) for s in boxes}
    return dict(frames=(avail[0], avail[-1]), n=nf, n_new=n_new,
                target=target, comp=comp, latest=latest,
                diff_latest=_latest_finite(dxd, dyd),
                png=os.path.join(out, PNG_NAME), npz=npz_path)


def _trailing_mean(arr, n=5):
    out = np.full_like(arr, np.nan, dtype=float)
    for i in range(len(arr)):
        seg = arr[max(0, i - n + 1):i + 1]
        seg = seg[np.isfinite(seg)]
        if len(seg) > 0:
            out[i] = np.mean(seg)
    return out


def _plot_timeseries(out, avail, dx, dy, dxd, dyd, obj_counts, sky_counts,
                     colors, diff_name, boxes, f_template, events_file,
                     last_n, show_all):
    nf = len(avail)
    AVG_N = 5
    vis = slice(nf - last_n, nf) if (not show_all and nf > last_n) else slice(None)
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
    oc_vis = []
    for star in boxes:
        col = colors[star]
        axes[0, 0].plot(avail, dx[star], "-o", ms=4, color=col, lw=0.9,
                        label=f"{star}  σ={np.nanstd(dx[star]):.3f}")
        axes[0, 0].plot(avail, _trailing_mean(dx[star], AVG_N), "-", color=col, lw=2, alpha=0.45)
        axes[1, 0].plot(avail, dy[star], "-o", ms=4, color=col, lw=0.9,
                        label=f"{star}  σ={np.nanstd(dy[star]):.3f}")
        axes[1, 0].plot(avail, _trailing_mean(dy[star], AVG_N), "-", color=col, lw=2, alpha=0.45)
        oc = obj_counts[star]
        med = float(np.nanmedian(oc[vis]))
        ocn = oc / med if med > 0 else oc
        oc_vis.append(ocn[vis])
        axes[2, 0].plot(avail, ocn, "-o", ms=4, color=col, lw=0.9, label=f"{star} med={med:.2e}")
        axes[2, 1].plot(avail, sky_counts[star], "-o", ms=4, color=col, lw=0.9,
                        label=f"{star} med={np.nanmedian(sky_counts[star]):.0f}")
    axes[0, 1].plot(avail, dxd, "-o", ms=4, color="k", lw=0.9, label=f"diff dx σ={np.nanstd(dxd):.3f}")
    axes[0, 1].plot(avail, _trailing_mean(dxd, AVG_N), "-", color="k", lw=2, alpha=0.45)
    axes[1, 1].plot(avail, dyd, "-o", ms=4, color="k", lw=0.9, label=f"diff dy σ={np.nanstd(dyd):.3f}")
    axes[1, 1].plot(avail, _trailing_mean(dyd, AVG_N), "-", color="k", lw=2, alpha=0.45)
    for ax in axes.flat:
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    for ax in (axes[0, 0], axes[1, 0], axes[0, 1], axes[1, 1]):
        ax.axhline(0, color="k", lw=0.4); ax.set_ylim(-1.5, 1.5)
    occ = np.concatenate(oc_vis) if oc_vis else np.array([])
    occ = occ[np.isfinite(occ)]
    y_lo_oc = max(0.0, min(0.98, float(np.nanmin(occ)) - 0.01)) if occ.size else 0.98
    axes[2, 0].axhline(1.0, color="k", lw=0.4); axes[2, 0].set_ylim(y_lo_oc, 1.02)
    axes[0, 0].set_title("dx (xcor on X profile) per star")
    axes[1, 0].set_title("dy (Gaussian fit to dip) per star")
    axes[0, 1].set_title(f"differential dx ({diff_name})")
    axes[1, 1].set_title(f"differential dy ({diff_name})")
    axes[2, 0].set_title("object counts (norm. to visible median)")
    axes[2, 1].set_title("sky pedestal (DN/pix)")
    axes[0, 0].set_ylabel("dx [px]"); axes[1, 0].set_ylabel("dy [px]")
    axes[2, 0].set_xlabel("frame"); axes[2, 1].set_xlabel("frame")
    events = load_events(events_file)
    if events:
        print(f"\nloaded {len(events)} motion events from {events_file}")
        for ax in axes.flat:
            mark_events(ax, events)
    if not show_all and nf > last_n:
        for ax in axes.flat:
            ax.set_xlim(avail[-last_n] - 0.5, avail[-1] + 0.5)
        tag = f"  [last {last_n} of {nf}]"
    else:
        tag = "  [all frames]"
    fig.suptitle(f"frames {avail[0]}..{avail[-1]} (n={nf})  template={f_template}{tag}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, PNG_NAME), dpi=120)
    plt.close(fig)
    print("wrote", os.path.join(out, PNG_NAME))


def _plot_box_montages(out, avail, boxes, box_imgs, dx, dy, mu_all, template_mu, f_template):
    nf = len(avail)
    for star, box in boxes.items():
        x_lo, x_hi = box["x_lo"], box["x_hi"]
        y_lo, y_hi = box["y_lo"], box["y_hi"]
        x_cen = 0.5 * (x_lo + x_hi - 1)
        y_cen_t = y_lo + template_mu[star]
        tpl = box_imgs[star][0]
        if tpl is None:
            continue
        vmin = float(np.nanpercentile(tpl, 5)); vmax = float(np.nanpercentile(tpl, 99))
        n_rows = (nf + N_COLS - 1) // N_COLS
        fig, axes = plt.subplots(n_rows, N_COLS, figsize=(2.0 * N_COLS, 2.4 * n_rows), squeeze=False)
        for i, fno in enumerate(avail):
            r, c = divmod(i, N_COLS); ax = axes[r, c]
            if box_imgs[star][i] is None:
                ax.axis("off"); continue
            ax.imshow(box_imgs[star][i], origin="lower",
                      extent=[x_lo - 0.5, x_hi - 0.5, y_lo - 0.5, y_hi - 0.5],
                      cmap="gray", vmin=vmin, vmax=vmax, interpolation="nearest", aspect="auto")
            cx = x_cen + dx[star][i]
            cy = y_lo + mu_all[star][i] if np.isfinite(mu_all[star][i]) else y_cen_t
            ax.axhline(cy, color="red", lw=0.7, alpha=0.85); ax.axvline(cx, color="red", lw=0.7, alpha=0.85)
            ax.set_title(f"f{fno} Δx={dx[star][i]:+.2f} Δy={dy[star][i]:+.2f}", fontsize=7)
            ax.set_xticks([]); ax.set_yticks([])
        for j in range(nf, n_rows * N_COLS):
            r, c = divmod(j, N_COLS); axes[r, c].axis("off")
        fig.suptitle(f"{star}  box {x_lo}..{x_hi} × {y_lo}..{y_hi} (template f{f_template})", fontsize=10)
        fig.tight_layout()
        p = os.path.join(out, f"box_montage_{cfgmod.sanitize(star)}.png")
        fig.savefig(p, dpi=120); plt.close(fig)
        print("wrote", p)
