"""
Shared core: the ~/guiding directory model, the observing-night date rule,
config.txt parsing (with safe fallback), and bundled-data access.

A "run dir" is one guided object on one night:
    ~/guiding/<night>-<target>-<comp>/
        config.txt          (settings)
        motion_events.txt   (keypress / nudge log)
        outputs/            (generated artifacts)

Nothing here is hardcoded to a user or absolute path: the data dir comes from
config.txt, the BPM ships inside the package, and the guiding root is
~/guiding (override with $HENRIETTA_GUIDING_DIR or --dir).
"""
from __future__ import annotations

import datetime as _dt
import os
import re
import shutil
from importlib import resources

# --- analysis defaults (used when config.txt is missing/unparseable) ---
DEFAULT_BOXES = {
    "TARGET": dict(x_lo=79, x_hi=122, y_lo=670, y_hi=852, role="target"),
    "COMP": dict(x_lo=879, x_hi=922, y_lo=670, y_hi=852, role="comp"),
}
DEFAULTS = dict(
    source="",
    bpm="bundled",
    start_frame=1,
    object="",
    filter="",
    pa=0.0,
    pixscale=0.163,   # guider plate scale, arcsec/pix
)


def bundled_bpm_path() -> str:
    """Filesystem path to the bad-pixel mask shipped inside the package."""
    return str(resources.files("henrietta_slit_guiding") / "data" / "bpm_25apr2026.fits")


def template_text() -> str:
    """The config.txt template shipped with the package."""
    return (resources.files("henrietta_slit_guiding") / "data" / "config_template.txt").read_text()


def guiding_root() -> str:
    """Root holding the per-object run dirs ($HENRIETTA_GUIDING_DIR or ~/guiding)."""
    return os.environ.get("HENRIETTA_GUIDING_DIR") or os.path.expanduser("~/guiding")


def night_date(now: _dt.datetime | None = None) -> str:
    """'YYYY_MM_DD' for the observing night.

    Observing-night rule: between 00:00 and 11:59 the night is the PREVIOUS
    calendar day (a post-midnight start still belongs to the night that began
    the evening before); from 12:00 on it is today.
    """
    now = now or _dt.datetime.now()
    d = now.date()
    if now.hour < 12:
        d = d - _dt.timedelta(days=1)
    return d.strftime("%Y_%m_%d")


def sanitize(name: str) -> str:
    """Filesystem-safe token for a directory-name component."""
    return re.sub(r"[^A-Za-z0-9._+-]", "_", name.strip())


def run_dir_name(target: str, comp: str, now: _dt.datetime | None = None) -> str:
    return f"{night_date(now)}-{sanitize(target)}-{sanitize(comp)}"


def roles_of(boxes: dict) -> tuple:
    """(target_names, comp_names) lists, by the role column."""
    targets = [n for n, b in boxes.items() if b.get("role", "").lower() == "target"]
    comps = [n for n, b in boxes.items() if b.get("role", "").lower() == "comp"]
    return targets, comps


def _default_config() -> dict:
    cfg = dict(DEFAULTS)
    cfg["bpm_path"] = bundled_bpm_path()
    cfg["boxes"] = {k: dict(v) for k, v in DEFAULT_BOXES.items()}
    return cfg


def load_config(path: str) -> dict:
    """Parse config.txt -> dict, with a safe fallback to built-in defaults.

    Keys returned: source, bpm (raw), bpm_path (resolved), start_frame,
    object, filter, pa, pixscale, boxes.  Recognised lines (a '#' starts a
    comment):
        source <dir>            bpm <bundled|path>
        start_frame <N>         object <substr>      filter <name>
        pa <deg>                pixscale <arcsec/pix>
        box <name> <role> <x_center> <x_halfwidth> <y_lo> <y_hi>
    Boxes need exactly one 'target' and one 'comp'.  ANY problem -> defaults
    + a printed warning, so the live loop never dies on a typo.
    """
    if not os.path.exists(path):
        return _default_config()
    cfg = dict(DEFAULTS)
    boxes = {}
    try:
        with open(path) as fh:
            for raw in fh:
                line = raw.split("#", 1)[0].strip()
                if not line:
                    continue
                kind, rest = (line.split(None, 1) + [""])[:2]
                kind = kind.lower()
                rest = rest.strip()
                if kind == "source":
                    cfg["source"] = rest
                elif kind == "bpm":
                    if rest:
                        cfg["bpm"] = rest
                elif kind == "start_frame":
                    cfg["start_frame"] = int(rest.split()[0])
                elif kind == "object":
                    cfg["object"] = rest
                elif kind == "filter":
                    cfg["filter"] = rest
                elif kind == "pa":
                    cfg["pa"] = float(rest.split()[0])
                elif kind == "pixscale":
                    cfg["pixscale"] = float(rest.split()[0])
                elif kind == "box":
                    parts = rest.split()
                    if len(parts) < 6:
                        raise ValueError(f"box needs 6 fields: {raw!r}")
                    name, role = parts[0], parts[1]
                    xc, hw, ylo, yhi = (int(round(float(p))) for p in parts[2:6])
                    if hw <= 0 or yhi <= ylo:
                        raise ValueError(f"bad box geometry: {raw!r}")
                    boxes[name] = dict(x_lo=xc - hw, x_hi=xc + hw + 1,
                                       y_lo=ylo, y_hi=yhi, role=role)
                else:
                    raise ValueError(f"unknown line (expected source/bpm/"
                                     f"start_frame/object/filter/pa/pixscale/"
                                     f"box): {raw!r}")
        targets, comps = roles_of(boxes)
        if len(targets) != 1 or len(comps) != 1:
            raise ValueError(f"need exactly one 'target' and one 'comp' role; "
                             f"got targets={targets} comps={comps}")
    except Exception as e:  # noqa: BLE001 - never let a bad config kill the loop
        print(f"  WARNING: config.txt unusable ({e}); using built-in defaults")
        return _default_config()
    # resolve the BPM: 'bundled' -> packaged file; else abs or relative-to-config
    if cfg["bpm"] == "bundled":
        cfg["bpm_path"] = bundled_bpm_path()
    elif os.path.isabs(cfg["bpm"]):
        cfg["bpm_path"] = cfg["bpm"]
    else:
        cfg["bpm_path"] = os.path.join(os.path.dirname(os.path.abspath(path)), cfg["bpm"])
    cfg["boxes"] = boxes
    return cfg


def target_comp(boxes: dict) -> tuple:
    """(target_name, comp_name) — assumes a validated config (one each)."""
    t, c = roles_of(boxes)
    return t[0], c[0]


def update_config_boxes(path: str, updates: dict) -> list:
    """Update (or append) `box` lines in config.txt for the given stars.

    `updates` maps name -> dict(role, x_center, halfwidth, y_lo, y_hi).  An
    existing `box <name> ...` line has its fields rewritten in place; a name
    with no line is appended.  A backup is written to `path + '.bak'` first.
    Returns [(name, 'updated'|'added'), ...].
    """
    def box_line(name, u):
        return (f"box  {name}  {u['role']}  {u['x_center']}  {u['halfwidth']}  "
                f"{u['y_lo']}  {u['y_hi']}   # x_center set by find-stars\n")

    with open(path) as fh:
        lines = fh.readlines()
    remaining = dict(updates)
    out, report = [], []
    for line in lines:
        parts = line.split("#", 1)[0].split()
        if len(parts) >= 2 and parts[0].lower() == "box" and parts[1] in remaining:
            name = parts[1]
            out.append(box_line(name, remaining.pop(name)))
            report.append((name, "updated"))
        else:
            out.append(line)
    for name, u in remaining.items():
        if out and not out[-1].endswith("\n"):
            out[-1] += "\n"
        out.append(box_line(name, u))
        report.append((name, "added"))

    shutil.copyfile(path, path + ".bak")
    with open(path, "w") as fh:
        fh.writelines(out)
    return report


def seeded_config_text(target: str, comp: str) -> str:
    """A config.txt pre-filled with the given names/roles (geometry/source
    left as EDIT placeholders).  Used by `init`."""
    return f"""\
# config.txt — {target} (target) + {comp} (comp).  Created by `init`.
# Fill in the EDIT fields, then `find-stars` for the box columns, then `watch`.
# Edit + save while watching and it redraws automatically.

source       /Users/henrietta/images/   # EDIT: directory holding hen####.fits
bpm          bundled             # 'bundled' = packaged mask, or a path
start_frame  1                   # EDIT: first frame number of the sequence
object       {target}            # substring matched in the OBJECT header
filter                           # EDIT: required FILTER (e.g. R-J); blank = any
pa           0.0                 # EDIT: instrument PA in deg (ROTANGLE) — keypress needs it
pixscale     0.163               # guider arcsec/pix, for the path/rate plots

# Extraction boxes — exactly one 'target' and one 'comp'.
#   box  <name>  <role>  <x_center>  <x_halfwidth>  <y_lo>  <y_hi>
# Set x_center with `find-stars`; keep y_lo..y_hi over the absorption feature
# the Y-fit tracks (R-J band: ~detector row 695).
box  {target}  target  100  21  670  852   # EDIT x_center via find-stars
box  {comp}  comp     900  21  670  852    # EDIT x_center via find-stars
"""
