# henrietta-slit-guiding — design

A `uv`-installable command-line tool for live slit-guiding analysis of a
two-star (target + comparison) spectroscopic field: it watches frames land,
measures each trace's cross-dispersion (X) and spectral (Y) drift relative
to a template, and tells the observer which TCS arrow keys re-center the
star — at the night's instrument PA.

This packages the tools developed in
`2026A/21June2026/HD136352-guiding/` into a portable command. **That live
directory is left untouched** (still mid-run); this is a separate build.

## Why

So a colleague (who uses `uv`) — or the next observer — can
`uv tool install` it on any Mac/Linux box and run it with no code editing:
all per-night / per-object settings live in a text `config.txt`, and the
bad-pixel mask ships inside the package.

## Distribution

- A real Python package `henrietta_slit_guiding` with a console script
  **`henrietta-slit-guiding`** (`[project.scripts]` in `pyproject.toml`).
- Install: `uv tool install <path-or-git-url>` (puts the command on PATH).
- Deps: `astropy numpy scipy matplotlib` (+ `argparse`/stdlib).
- Bundled package data: the BPM (`bpm_25apr2026.fits`) and a
  `config_template.txt`. Accessed via `importlib.resources` — no absolute
  paths anywhere.

## Runtime data model

Each guided object on a night is its own self-contained directory under
`~/guiding/`:

```
~/guiding/
  2026_06_21-HD136352-nu01Lup/      <- one object, one night
    config.txt                      <- settings (you edit)
    motion_events.txt               <- keypress / nudge log (you edit)
    outputs/                        <- generated: dx_dy png, live.html,
                                       box_overlay.png, .npz, motion_cache.json
  2026_06_21-WASP123-HD200000/      <- another object, same night
    ...
```

- **Dir name** = `<night>-<target>-<comp>`; `<target>`/`<comp>` come from the
  `box … target` / `box … comp` lines in the config (names are free; the
  code keys on the `role`). *(Decision: target-comp, not target-only.)*
- **`<night>`** = local calendar date with the observing-night rule:
  **00:00–11:59 → previous day**, otherwise today. (So a post-midnight start
  still belongs to the night that began the evening before.)
- `~/guiding` is the default root; override with `--dir <path>` or
  `$HENRIETTA_GUIDING_DIR`.

## config.txt schema

Plain text; `#` comments; blank lines ignored. `init <target> <comp>` writes
this file **pre-seeded** — the box names/roles and `object` come from the
arguments; `source` and the box columns are placeholders you fill in. Bad /
half-saved file → the command falls back to built-in defaults and warns
(never crashes the loop).

```
source       /Volumes/Extreme Pro/HenJune2026   # dir of hen####.fits (spaces ok)
bpm          bundled                             # 'bundled' = packaged BPM, or a path
start_frame  83                                  # first frame of the sequence
object       136352                              # substring required in OBJECT header
filter       R-J                                 # required FILTER (blank = any)
# box  <name>  <role>  <x_center>  <x_halfwidth>  <y_lo>  <y_hi>
#   exactly one 'target' + one 'comp'; R-J absorption feature ~row 695,
#   so keep y_lo ~640-670.  `find-stars` prints starter values.
box  HD136352  target  109   21  670  852
box  nu01Lup   comp    1943  21  670  852
```

## Subcommands

All commands operate on the **config.txt in the current directory** (so you
`cd` into an object's dir to act on it); `--dir` overrides.

- **`init <target> <comp>`** — single pass. The **directory comes first and
  seeds the config**: it computes the night date, creates
  `~/guiding/<night>-<target>-<comp>/`, and writes a `config.txt` there
  **auto-populated from the arguments** — the two `box` lines already carry
  the target/comp names + roles (with placeholder geometry to fill in), and
  `object` is seeded with the target name. Adds an empty `motion_events.txt`
  and prints the `cd`. You then fill in `source` and the box columns (run
  `find-stars` to get the columns). Names with spaces are sanitized for the
  dir name but kept verbatim in the config. *(Decision: single-pass, dir
  drives config — not "config drives dir".)*
- **`watch`** — the live loop (explicit verb, not a bare default). Regenerates
  `outputs/` on each new frame / `motion_events.txt` edit / `config.txt`
  edit; writes `box_overlay.png` once at start + on config change; serves the
  auto-refreshing `live.html`. **Prints a live status readout each redraw**
  (frames, per-star + differential drift, recommended keys) so the terminal
  is useful over ssh without the browser. *(Decision: verb = `watch`.)*
- **`keypress [--dx --dy] [--pa --scale]`** — drift → TCS arrow presses;
  default reads the latest drift for the target from `outputs/*.npz`.
- **`find-stars`** — Y-smash the first frame, peak-find the trace columns,
  print suggested `box` lines for `config.txt`.
- **`headers`** — dump OBJECT/FILTER/AIRMASS/ROTANGLE/DATE-OBS for frames in
  range, to pick `start_frame` and confirm the `object` string.
- **`overlay`** — (re)write `outputs/box_overlay.png`.
- **`path`** / **`rate`** — post-night plots (guider keypress trail; commanded
  motion histogram), reading `motion_events.txt`.

### `watch` status readout (the "mini-TUI", no deps)

After each redraw, print a compact block — keeps 80% of a TUI's value for ~5%
of the cost; a real Textual TUI is an explicit **v2** if logging-by-editor
proves annoying in use:

```
21:42:10  frames 83..208 (n=126)   new:1 cached:125
  HD136352  dx -0.08  dy -0.46        nu01Lup  dx +0.01  dy +0.00
  diff      dx -0.09  dy -0.46   |d|=0.47 px
  ==> PRESS:  UP 1   LEFT 1
```

## Carried-over behavior (already built, just repackaged)

- **Incremental fit cache** (`outputs/motion_cache.json`): only new frames are
  fit; reuse keyed on file-size + signature (boxes/template/params). `--no-cache`.
- **`live.html`**: browser page that re-fetches the PNG every 2 s (Preview
  won't auto-refresh a background window); has a "box layout" link.
- **Names fully from config** (role-driven); `.npz` uses dynamic `dx_<name>`
  keys + `target`/`comp` metadata.
- **Self-locating paths** everywhere; nothing hardcoded to a user/home.

## Package layout

```
guider-tools/                         (repo; install from here)
  pyproject.toml                      ([project.scripts], deps, package-data)
  README.md  DESIGN.md
  src/henrietta_slit_guiding/
    __init__.py
    cli.py            argparse subcommands -> dispatch
    config.py         load_config, night-date, run-dir resolution, paths
    montage.py        analysis + plot + cache + box overlay (was motion_montage2)
    watch.py          watch loop + status readout
    keypress.py  findstars.py  headers.py  plots.py   (path/rate)
    data/  bpm_25apr2026.fits   config_template.txt
  tests/  test_config.py  test_watch.py  test_montage_cache.py  ...
```

## Testing

- Pure-function unit tests (no data): config parse + fallbacks, night-date
  rule (incl. the noon boundary), run-dir name derivation, cache
  round-trip/invalidation, watch change-detection, keypress matrix math.
- Integration smoke against the live drive (it's mounted): `init` → edit →
  `watch` one pass → assert `outputs/` artifacts + non-NaN drift; `keypress`
  reads the npz; cached == `--no-cache` (byte-identical).

## Out of scope (for now)

- Full Textual **TUI** — deferred to v2 (revisit once the CLI is in use).
- Auto-recenter / box tracking — explicitly not wanted; box stays fixed.
- Touching the live `2026A/21June2026/HD136352-guiding/` directory.
- Generalizing beyond two stars (one target + one comp).
```
