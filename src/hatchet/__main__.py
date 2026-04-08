# src/hatchet/__main__.py
import argparse
from importlib.metadata import version

from hatchet.cluster_bins.cluster_bins import run as hatchet_cluster_bins
from hatchet.compute_cn.compute_cn import run as hatchet_compute_cn
from hatchet.plot.plot_cn import run as hatchet_plot_cn
from hatchet.plot.plot_cnp_panel import run as hatchet_plot_panel
from hatchet.evaluate.evaluate import run as hatchet_evaluate
from hatchet.hatchet_parser import *
from hatchet.utils import setup_logging, log_arguments


def main(argv=None):
    parser = argparse.ArgumentParser(prog="hatchet")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {version('hatchet')}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_clu = subparsers.add_parser("cluster-bins", help="run cluster-bins")
    add_arguments_cluster_bins(p_clu)
    p_clu.set_defaults(func=hatchet_cluster_bins)

    p_ccn = subparsers.add_parser("compute-cn", help="run hatchet_compute_cn")
    add_arguments_compute_cn(p_ccn)
    p_ccn.set_defaults(func=hatchet_compute_cn)

    p_plot_cn = subparsers.add_parser("plot-cn", help="run plot-cn")
    add_arguments_plot_cn(p_plot_cn)
    p_plot_cn.set_defaults(func=hatchet_plot_cn)

    p_plot_cnp = subparsers.add_parser("plot-panel", help="run plot-panel")
    add_arguments_plot_panel(p_plot_cnp)
    p_plot_cnp.set_defaults(func=hatchet_plot_panel)

    p_eval = subparsers.add_parser(
        "evaluate", help="evaluate CN solutions against somatic SNVs"
    )
    add_arguments_evaluate(p_eval)
    p_eval.set_defaults(func=hatchet_evaluate)

    args = parser.parse_args(argv)
    setup_logging(args)
    log_arguments(args)
    args.func(args)


if __name__ == "__main__":
    main()
