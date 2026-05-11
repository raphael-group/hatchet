# src/hatchet/__main__.py
import argparse
from importlib.metadata import version

from hatchet.cluster_bins.cluster_bins import run as hatchet_cluster_bins
from hatchet.compute_cn.compute_cn import run as hatchet_compute_cn
from hatchet.plot.plot_cn import run as hatchet_plot_cn
from hatchet.plot.plot_cnp_panel import run as hatchet_plot_panel
from hatchet.evaluate.evaluate import run as hatchet_evaluate
from hatchet.hatchet_parser import (
    add_arguments_cluster_bins,
    add_arguments_compute_cn,
    add_arguments_evaluate,
    add_arguments_plot_cn,
    add_arguments_plot_panel,
)


SUBCOMMANDS = [
    ("cluster-bins", add_arguments_cluster_bins, hatchet_cluster_bins),
    ("compute-cn", add_arguments_compute_cn, hatchet_compute_cn),
    ("plot-cn", add_arguments_plot_cn, hatchet_plot_cn),
    ("plot-panel", add_arguments_plot_panel, hatchet_plot_panel),
    ("evaluate", add_arguments_evaluate, hatchet_evaluate),
]


def main(argv=None):
    parser = argparse.ArgumentParser(prog="hatchet")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {version('hatchet')}"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    for name, add_args, fn in SUBCOMMANDS:
        p = sub.add_parser(name, help=f"run {name}")
        add_args(p)
        p.set_defaults(func=fn)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
