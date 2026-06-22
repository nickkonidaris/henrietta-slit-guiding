# henrietta-slit-guiding

Live slit-guiding analysis for a two-star (target + comparison) spectroscopic
field on the Henrietta / Swope instrument. It watches frames land, measures
each trace's cross-dispersion (X) and spectral (Y) drift relative to a
template, and tells you which TCS arrow keys re-center the star at the night's
instrument PA.

Install once, run anywhere — no copying code per machine or per night. See
[`DESIGN.md`](DESIGN.md) for the full design rationale.

## Install (uv)

```sh
uv tool install git+https://github.com/nickkonidaris/henrietta-slit-guiding
# later:  uv tool upgrade henrietta-slit-guiding
```
(Plain pip works too: `pipx install git+...` or `pip install git+...`.)

The bad-pixel mask ships inside the package — nothing else to fetch.

## Workflow

Each guided object on a night is its own directory under `~/guiding/`,
holding its `config.txt`, keypress log, and outputs.

```sh
henrietta-slit-guiding init HD136352 nu01Lup
#  -> ~/guiding/2026_06_21-HD136352-nu01Lup/   (config.txt seeded with the names)

cd ~/guiding/2026_06_21-HD136352-nu01Lup
#  edit config.txt: source, start_frame, filter, pa
henrietta-slit-guiding headers       # confirm OBJECT/FILTER, pick start_frame
henrietta-slit-guiding find-stars    # finds the columns and writes them into config.txt
henrietta-slit-guiding watch         # the live loop (auto-opens live.html in your browser)
```

`watch` opens `outputs/live.html` for you (a browser page that auto-refreshes
the plot — macOS Preview won't reload a background window); pass `--no-open` to
skip that. `find-stars` writes the box columns straight into `config.txt`
(backing it up to `config.txt.bak`), then redraws and **opens
`box_overlay.png`** so you can check the target/comp assignment on the actual
frame (`--no-open` to skip, `--dry-run` to only print). Log nudges in
`motion_events.txt`; the plot redraws on each new frame / events edit / config
edit. `watch` also prints a status block each redraw:

```
21:42:10  frames 83..208 (n=126)   new:1 cached:125
  HD136352     dx -0.08  dy -0.46        nu01Lup      dx +0.01  dy +0.00
  diff         dx -0.09  dy -0.46   |d|=0.47 px
  ==> PRESS:  UP 1   LEFT 1
```

### Other subcommands

```sh
henrietta-slit-guiding keypress              # drift -> arrow keys (or --dx .. --dy ..)
henrietta-slit-guiding find-stars            # write columns to config.txt + open box_overlay.png (--dry-run / --no-open)
henrietta-slit-guiding overlay               # (re)write outputs/box_overlay.png
henrietta-slit-guiding path                  # keypress-trail plot
henrietta-slit-guiding rate                  # commanded-motion histogram
```

All commands act on the `config.txt` in the current directory (or `--dir`).
The guiding root is `~/guiding` (override with `$HENRIETTA_GUIDING_DIR`).

## config.txt

```
source       /Users/henrietta/images/           # dir of hen####.fits
bpm          bundled                             # packaged mask, or a path
start_frame  83
object       136352                              # substring required in OBJECT
filter       R-J                                 # blank = any
pa           9.336                               # instrument PA (ROTANGLE) deg
pixscale     0.163                               # guider arcsec/pix
# box <name> <role> <x_center> <x_halfwidth> <y_lo> <y_hi>  (one target, one comp)
box HD136352 target 109  21 670 852
box nu01Lup  comp    1943 21 670 852
```

Names are yours; the code keys on the `role`. The R-J absorption feature the
Y-fit tracks sits near detector row ~695, so keep `y_lo`~640–670. A bad or
half-saved line is ignored (built-in defaults used) so editing during `watch`
can't crash the loop.

## Tests

```sh
python -m pytest            # or: python tests/test_config.py
```
