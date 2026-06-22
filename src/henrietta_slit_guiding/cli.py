"""
henrietta-slit-guiding — command-line entry point.

Subcommands:
  init <target> <comp>   create ~/guiding/<night>-<target>-<comp>/ with a
                         seeded config.txt + empty motion_events.txt
  watch                  the live loop (status readout + auto-refresh plot)
  keypress               drift -> TCS arrow presses
  find-stars             suggest box columns for config.txt
  headers                dump OBJECT/FILTER/ROTANGLE to pick start_frame
  overlay                (re)write outputs/box_overlay.png
  path / rate            post-night plots from motion_events.txt

All commands except `init` act on the config.txt in the current directory
(or --dir).
"""
from __future__ import annotations

import argparse
import os
import sys

from . import config as cfgmod


def _resolve_run_dir(args) -> str:
    return os.path.abspath(args.dir) if args.dir else os.getcwd()


def _load_or_die(run_dir: str) -> dict:
    cfg_path = os.path.join(run_dir, "config.txt")
    if not os.path.exists(cfg_path):
        sys.exit(f"no config.txt in {run_dir}\n"
                 f"  cd into a run dir (made by `init`), or pass --dir.")
    return cfgmod.load_config(cfg_path)


def _cmd_init(args):
    target, comp = args.target, args.comp
    root = cfgmod.guiding_root()
    name = cfgmod.run_dir_name(target, comp)
    run_dir = os.path.join(root, name)
    cfg_path = os.path.join(run_dir, "config.txt")
    if os.path.exists(cfg_path):
        print(f"already exists: {cfg_path}\n(leaving it untouched)")
    else:
        os.makedirs(os.path.join(run_dir, "outputs"), exist_ok=True)
        with open(cfg_path, "w") as fh:
            fh.write(cfgmod.seeded_config_text(target, comp))
        events = os.path.join(run_dir, "motion_events.txt")
        if not os.path.exists(events):
            with open(events, "w") as fh:
                fh.write("# <frame> <direction up/down/left/right> <amount> "
                         "[comment]\n")
    print(f"created {run_dir}")
    print(f"\nnext:")
    print(f"  cd {run_dir}")
    print(f"  # edit config.txt: source, start_frame, filter, pa")
    print(f"  henrietta-slit-guiding find-stars   # then paste box columns")
    print(f"  henrietta-slit-guiding watch")


def _cmd_watch(args):
    from . import watch
    run_dir = _resolve_run_dir(args)
    _load_or_die(run_dir)  # validate config exists/parses now
    watch.watch(run_dir, last_n=args.last_n, interval=args.interval,
                open_browser=not args.no_open)


def _cmd_keypress(args):
    from . import keypress
    run_dir = _resolve_run_dir(args)
    cfg = _load_or_die(run_dir)
    keypress.run_keypress(run_dir, cfg, dx=args.dx, dy=args.dy, scale=args.scale,
                          deadband=args.deadband, star=args.star, pa=args.pa)


def _cmd_find_stars(args):
    from . import findstars
    run_dir = _resolve_run_dir(args)
    findstars.find_stars(_load_or_die(run_dir), run_dir, dry_run=args.dry_run)


def _cmd_headers(args):
    from . import headers
    headers.dump_headers(_load_or_die(_resolve_run_dir(args)))


def _cmd_overlay(args):
    from . import montage
    run_dir = _resolve_run_dir(args)
    montage.write_box_overlay(run_dir, _load_or_die(run_dir))


def _cmd_path(args):
    from . import plots
    run_dir = _resolve_run_dir(args)
    plots.plot_path(run_dir, _load_or_die(run_dir))


def _cmd_rate(args):
    from . import plots
    run_dir = _resolve_run_dir(args)
    plots.plot_rate(run_dir, _load_or_die(run_dir))


def build_parser():
    p = argparse.ArgumentParser(
        prog="henrietta-slit-guiding",
        description="Live slit-guiding analysis for a two-star field.")
    p.add_argument("--dir", default=None,
                   help="Run directory (default: current directory).")
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("init", help="create a new run dir from star names")
    pi.add_argument("target")
    pi.add_argument("comp")
    pi.set_defaults(func=_cmd_init)

    pw = sub.add_parser("watch", help="live loop")
    pw.add_argument("--last-n", type=int, default=30)
    pw.add_argument("--interval", type=float, default=5.0)
    pw.add_argument("--no-open", action="store_true",
                    help="don't auto-open live.html in the browser")
    pw.set_defaults(func=_cmd_watch)

    pk = sub.add_parser("keypress", help="drift -> TCS arrow presses")
    pk.add_argument("--dx", type=float, default=None)
    pk.add_argument("--dy", type=float, default=None)
    pk.add_argument("--scale", type=float, default=1.0)
    pk.add_argument("--deadband", type=float, default=0.2)
    pk.add_argument("--star", default=None, help="default: the config target")
    pk.add_argument("--pa", type=float, default=None, help="default: config pa")
    pk.set_defaults(func=_cmd_keypress)

    pf = sub.add_parser("find-stars", help="find box columns and write them to config.txt")
    pf.add_argument("--dry-run", action="store_true",
                    help="print only; don't modify config.txt")
    pf.set_defaults(func=_cmd_find_stars)
    sub.add_parser("headers", help="dump frame headers").set_defaults(func=_cmd_headers)
    sub.add_parser("overlay", help="write box_overlay.png").set_defaults(func=_cmd_overlay)
    sub.add_parser("path", help="keypress-trail plot").set_defaults(func=_cmd_path)
    sub.add_parser("rate", help="commanded-motion histogram").set_defaults(func=_cmd_rate)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
