import logging
import argparse


def solver_available(solver_type: str):
    from pyomo import environ as pe

    found_solver = False
    if solver_type == "gurobi":
        found_solver = pe.SolverFactory("gurobi", solver_io="python").available(
            exception_flag=False
        )
    else:
        found_solver = pe.SolverFactory(solver_type).available(exception_flag=False)
    if found_solver:
        logging.info(f"solver={solver_type} is available.")
    return found_solver


##################################################
def add_arguments_cluster_bins(parser: argparse.ArgumentParser):
    ##################################################
    # general inputs for bulk data
    parser.add_argument(
        "--bb_dir",
        required=True,
        type=str,
        help="Directory containing NPZ count matrices and bb.tsv.gz input files (default: bb)",
    )
    parser.add_argument(
        "--bbc_dir",
        required=True,
        type=str,
        help="Output directory for BBC/seg TSV files",
    )

    ##################################################
    # Parameters
    parser.add_argument(
        "--minK",
        required=False,
        default=3,
        type=int,
        help="Minimum number of HMM cluster states (default: 3)",
    )
    parser.add_argument(
        "--maxK",
        required=False,
        default=30,
        type=int,
        help="Maximum number of HMM cluster states (default: 30)",
    )
    parser.add_argument(
        "-t",
        required=False,
        default=1e-6,
        type=float,
        help="initial off-diagonal transition mass (default: 1e-6)",
    )
    parser.add_argument(
        "--restarts",
        required=False,
        default=30,
        type=int,
        help="#restarts per K (default: 30)",
    )
    parser.add_argument(
        "--top_restarts",
        required=False,
        default=30,
        type=int,
        help="number of top-scoring inits to run full EM on (default: 30)",
    )
    parser.add_argument(
        "--n_local_trials",
        required=False,
        default=3,
        type=int,
        help="Number of candidate bins to evaluate per k-means++ seeding step; "
        "best by log-likelihood is kept (default: 3, 1 = original single-draw)",
    )
    parser.add_argument(
        "--niters",
        required=False,
        default=50,
        type=int,
        help="Number of EM iterations per restart (default: 50)",
    )

    parser.add_argument(
        "--decode_method",
        required=False,
        choices=["viterbi", "map"],
        type=str,
        help="HMM decoding method: viterbi (most-likely path) or map (marginal per-bin posterior) (default: map)",
        default="map",
    )
    parser.add_argument(
        "--score_method",
        required=False,
        choices=["bic", "icl"],
        type=str,
        help="Model selection criterion (default: icl)",
        default="icl",
    )
    parser.add_argument(
        "--init_method",
        required=False,
        choices=["cna_plus_plus", "kmeans_plus_plus"],
        type=str,
        default="cna_plus_plus",
        help="HMM initialization method: 'cna_plus_plus' (HMM-aware seeding) or 'kmeans_plus_plus' (sklearn KMeans++ in [RDR, mhBAF] space) (default: cna_plus_plus)",
    )

    parser.add_argument(
        "--bb_quantile",
        required=False,
        default=0.2,
        type=float,
        help="Fraction of near-balanced bins (BAF ≈ 0.5) used to initialize Beta-Binomial dispersion (default: 0.2)",
    )

    parser.add_argument(
        "--min_tau",
        required=False,
        default=50,
        type=float,
        help="Minimum Beta-Binomial dispersion tau (default: 50)",
    )

    parser.add_argument(
        "--max_tau",
        required=False,
        default=500,
        type=float,
        help="Maximum Beta-Binomial dispersion tau (default: 500)",
    )

    parser.add_argument(
        "--baf_eps",
        required=False,
        default=1e-3,
        type=float,
        help="BAF mean Brent search bounds [baf_eps, 1-baf_eps]; related to sequencing error floor (default: 1e-3)",
    )

    parser.add_argument(
        "--min_covar",
        required=False,
        default=1e-3,
        type=float,
        help="Minimum RDR variance floor applied after each M-step (default: 1e-3)",
    )

    parser.add_argument(
        "--ig_alpha",
        required=False,
        default=10.0,
        type=float,
        help="Inverse-gamma prior shape parameter for RDR variance updates (default: 10.0)",
    )

    parser.add_argument(
        "--tau_iters",
        required=False,
        default=0,
        type=int,
        help="Number of EM iterations during which BAF dispersion tau is updated (default: 0)",
    )

    parser.add_argument(
        "--seed",
        required=False,
        default=42,
        type=int,
        help="random seed for HMM init step (default: 42)",
    )

    parser.add_argument(
        "--log_rdr",
        action="store_true",
        default=False,
        help="Use log(RDR) instead of raw RDR in the Gaussian emission (default: False).",
    )

    ##################################################
    # aux files
    parser.add_argument(
        "--genome_size",
        required=True,
        type=str,
        help="Reference chromosome sizes file (e.g., hg19.chrom.sizes)",
    )
    parser.add_argument(
        "--verbosity",
        required=False,
        default=0,
        type=int,
        help="verbose level, 0, 1, 2 (default: 0)",
    )
    return parser


##################################################
def add_arguments_compute_cn(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--result_dir",
        required=True,
        type=str,
        help="Output directory for computed CN results (default: results)",
    )
    parser.add_argument(
        "--bbc",
        required=True,
        type=str,
        help="Filename for BBC table (e.g., results/best.bbc.ucn)",
    )
    parser.add_argument(
        "--seg",
        required=True,
        type=str,
        help="Filename for SEG table (e.g., results/best.seg.ucn)",
    )

    ##################################################
    parser.add_argument(
        "--mode",
        required=False,
        choices=["both", "cd", "ilp"],
        type=str,
        help="Solver mode (default: ilp)",
        default="ilp",
    )

    parser.add_argument(
        "--solver",
        required=False,
        choices=["gurobi", "cbc"],
        type=str,
        help="ILP solver (default: gurobi)",
        default="gurobi",
    )

    parser.add_argument(
        "--timelimit",
        required=False,
        default=None,
        type=int,
        help="ILP solver timelimit in seconds (default: None)",
    )

    parser.add_argument(
        "--pool_size",
        required=False,
        default=1,
        type=int,
        help="Number of Gurobi solution-pool solutions to collect (mode 0, default: 1 = disabled). "
        "Only effective with --solver gurobi.",
    )

    parser.add_argument(
        "--pool_gap",
        required=False,
        default=None,
        type=float,
        help="Relative optimality gap for Gurobi solution pool (default: None = keep all). "
        "E.g. 0.0 keeps only optimal, 0.1 keeps within 10%% of optimal.",
    )
    ##################################################
    # preprocessing
    parser.add_argument(
        "--filter_cluster",
        action="store_true",
        default=False,
        required=False,
        help="Enable RDR/BAF variance outlier filtering of clusters before optimization (default: false)",
    )
    parser.add_argument(
        "--filter_std",
        required=False,
        default=2.0,
        type=float,
        help="Filter clusters whose variance deviates from mean by filter_std * std(variances) (default: 2.0)",
    )
    parser.add_argument(
        "--min_nbins",
        required=False,
        default=10,
        type=int,
        help="Minimum number of bins a cluster must have to be retained (default: 10)",
    )
    parser.add_argument(
        "--ub_nbins",
        required=False,
        default=50,
        type=int,
        help="Upper bound on bins: variance-outlier filtering only applies to clusters with #bins <= ub_nbins (default: 50)",
    )

    parser.add_argument(
        "--segment",
        action="store_true",
        default=False,
        required=False,
        help="Use genomic-segment-level data instead of cluster-level summaries (default: false)",
    )

    parser.add_argument(
        "--bal_tost_margin",
        required=False,
        default=3e-2,
        type=float,
        help="BAF-tolerance, locate balanced clusters (default: 0.03)",
    )

    parser.add_argument(
        "--bal_tost_alpha",
        type=float,
        required=False,
        default=0.05,
        help="TOST equivalence test significance level (default: 0.05)",
    )
    parser.add_argument(
        "--tol_nstd",
        type=float,
        required=False,
        default=1.0,
        help="Number of std deviations for CN scoring tolerance (default: 1.0)",
    )
    parser.add_argument(
        "--fcn_ci_alpha",
        type=float,
        required=False,
        default=0.05,
        help="Significance level for fractional CN confidence intervals (default: 0.05)",
    )

    ##################################################
    # model parameters
    parser.add_argument(
        "--minClone",
        required=False,
        default=2,
        type=int,
        help="Minimum number of tumor clones to solve for (default: 2)",
    )
    parser.add_argument(
        "--maxClone",
        required=False,
        default=4,
        type=int,
        help="Maximum number of tumor clones to solve for (default: 4)",
    )

    parser.add_argument(
        "--diploid",
        action="store_true",
        default=False,
        required=False,
        help="Solve under diploid assumption (cn_max=6)",
    )
    parser.add_argument(
        "--tetraploid",
        action="store_true",
        default=False,
        required=False,
        help="Solve under tetraploid/WGD assumption (cn_max=12)",
    )

    ##################################################
    # constraints & regularization
    parser.add_argument(
        "--reg_term",
        required=False,
        choices=["RAW", "MAXCN", "DROOT_SUM", "DADJ_SUM"],
        type=str,
        help="regularization term (default: MAXCN)",
        default="MAXCN",
    )
    parser.add_argument(
        "--reg_steps",
        required=False,
        default=10,
        type=int,
        help="Number of steps in the regularization path (default: 10)",
    )
    parser.add_argument(
        "--reg_stepsize",
        required=False,
        default=0.01,
        type=float,
        help="Multiplicative step size between regularization path values (default: 0.01)",
    )

    parser.add_argument(
        "--no_ampdel",
        action="store_true",
        default=False,
        required=False,
        help="Disable the amp/del symmetry constraint (default: off)",
    )
    parser.add_argument(
        "--num_cnstates",
        required=False,
        default=-1,
        type=int,
        help="Constrain the number of distinct CN states per clone (-1 = unconstrained, default: -1)",
    )
    parser.add_argument(
        "-eD",
        "--diploidcmax",
        type=int,
        required=False,
        default=6,
        help=(
            "Maximum copy-number value overall segments (default: 6, 0 means inferred from scaled fractional copy "
            "numbers)"
        ),
    )
    parser.add_argument(
        "-eT",
        "--tetraploidcmax",
        type=int,
        required=False,
        default=12,
        help=(
            "Maximum copy-number value overall segments (default: 12, 0 means inferred from scaled fractional "
            "copy numbers)"
        ),
    )
    parser.add_argument(
        "--min_prop",
        required=False,
        default=0.01,
        type=float,
        help="minimum clone proportion (default: 0.01)",
    )

    parser.add_argument(
        "--purities",
        required=False,
        default=None,
        type=lambda s: {
            k: float(v) for k, v in (pair.split(":") for pair in s.split(";"))
        },
        help=(
            "Semicolon-separated sample:purity pairs (e.g. sample1:0.80;sample2:0.70). "
            "When provided, each sample's normal-clone proportion is fixed to 1 - purity in the ILP solver."
        ),
    )

    ##################################################
    # coordinate descent parameters
    parser.add_argument(
        "--cd_niters",
        required=False,
        default=10,
        type=int,
        help="CD: max outer CD iterations per seed (default: 10)",
    )
    parser.add_argument(
        "--cd_convergence_iters",
        required=False,
        default=2,
        type=int,
        help="CD: consecutive convergence iterations required to stop (default: 2)",
    )
    parser.add_argument(
        "--cd_nseeds",
        required=False,
        default=400,
        type=int,
        help="CD: number of random restarts (default: 400)",
    )
    parser.add_argument(
        "--cd_njobs",
        required=False,
        default=8,
        type=int,
        help="CD: number of parallel worker processes (default: 8)",
    )
    parser.add_argument(
        "--cd_seed",
        required=False,
        default=42,
        type=int,
        help="CD: random seed for reproducibility (default: 42)",
    )

    parser.add_argument(
        "--verbosity",
        required=False,
        default=0,
        type=int,
        help="verbose level, 0, 1, 2 (default: 0)",
    )
    parser.add_argument(
        "--genome_size",
        required=True,
        type=str,
        help="Reference chromosome sizes file",
    )
    parser.add_argument(
        "--region_bed",
        required=True,
        type=str,
        help="Reference chromosome BED file",
    )
    return parser


##################################################
def parse_arguments_compute_cn(argv=None):
    parser = argparse.ArgumentParser(description="HATCHet compute-cn")
    add_arguments_compute_cn(parser)
    args = parser.parse_args(argv)

    # Validate solver availability only when ILP is needed
    if args.mode in ("ilp", "both"):
        if not solver_available(args.solver):
            if args.solver == "gurobi":
                parser.error(
                    "Gurobi solver is not available. "
                    "Ensure gurobipy is installed and a valid Gurobi license is active "
                    "(check GRB_LICENSE_FILE or ~/.gurobi/gurobi.lic)."
                )
            else:
                parser.error(
                    f"Solver '{args.solver}' is not available. "
                    "Ensure the corresponding Pyomo solver backend is installed and on PATH."
                )
    return args


##################################################
def add_arguments_plot_cn(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--bbc",
        required=True,
        type=str,
        help="Filename for BBC table (e.g., results/best.bbc.ucn)",
    )
    parser.add_argument(
        "--seg",
        required=True,
        type=str,
        help="Filename for SEG table (e.g., results/best.seg.ucn)",
    )
    parser.add_argument(
        "--genome_size",
        required=True,
        type=str,
        help="Reference chromosome sizes file (e.g., hg19.chrom.sizes)",
    )
    parser.add_argument(
        "--region_bed",
        required=True,
        type=str,
        help="Reference chromosome BED file (e.g., hg19.chrom.bed)",
    )
    parser.add_argument(
        "-g",
        "--gamma_file",
        required=True,
        type=str,
        help="gamma scaling factor file, output of HATCHet compute-cn",
    )
    parser.add_argument(
        "-s",
        "--solfile",
        required=False,
        type=str,
        help="Optional solution file to override CN states in BBC table",
        default=None,
    )
    parser.add_argument(
        "-O",
        "--plot_dir",
        required=True,
        type=str,
        help="Directory for output files",
    )
    parser.add_argument(
        "--dpi",
        required=False,
        type=int,
        help="image resolution (default: 300)",
        default=300,
    )
    parser.add_argument(
        "--img_type",
        required=False,
        choices=["pdf", "png", "svg"],
        type=str,
        help="file format (default: png)",
        default="png",
    )
    parser.add_argument(
        "--transparent",
        required=False,
        action="store_true",
        default=False,
        help="transparent background (default: False)",
    )
    parser.add_argument(
        "--keep_gap",
        required=False,
        action="store_true",
        default=False,
        help="keep gap region in the plot (default: False)",
    )
    parser.add_argument(
        "--tail_alpha",
        required=False,
        type=float,
        help="transparency on the tail region per CN state (default: 0.8)",
        default=0.8,
    )
    parser.add_argument(
        "--center_alpha",
        required=False,
        type=float,
        help="transparency on the center region per CN state (default: 1.0)",
        default=1.0,
    )
    parser.add_argument(
        "--onetail_area",
        required=False,
        type=float,
        help="area for each tail per CN state to set transparency (default: 0.025)",
        default=0.025,
    )
    parser.add_argument(
        "--maxlim_fcn",
        required=False,
        type=int,
        help="figure axis limit for FCN (default: 30)",
        default=30,
    )
    parser.add_argument(
        "--ploidy",
        required=True,
        choices=["diploid", "tetraploid"],
        help="Ploidy of the solution (selects gamma column from gamma file).",
    )
    return parser


##################################################
def add_arguments_plot_panel(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--panel_file",
        required=True,
        type=str,
        help="Panel TSV file listing per-sample BBC UCN paths",
    )
    parser.add_argument(
        "--genome_size",
        required=True,
        type=str,
        help="Reference chromosome sizes file (e.g., hg19.chrom.sizes)",
    )
    parser.add_argument(
        "--region_bed",
        required=True,
        type=str,
        help="Reference chromosome BED file (e.g., hg19.chrom.bed)",
    )
    parser.add_argument(
        "--width",
        required=False,
        type=int,
        help="panel image width (default: 20)",
        default=20,
    )
    parser.add_argument(
        "--height",
        required=False,
        type=int,
        help="panel image height per row (default: 1)",
        default=1,
    )
    parser.add_argument(
        "--show_clone_name",
        required=False,
        action="store_true",
        default=False,
        help="plot clone name (default: False)",
    )
    parser.add_argument(
        "--show_prop",
        required=False,
        action="store_true",
        default=False,
        help="plot clone proportion (default: False)",
    )
    parser.add_argument(
        "--dpi",
        required=False,
        type=int,
        help="image resolution (default: 300)",
        default=300,
    )
    parser.add_argument(
        "--transparent",
        required=False,
        action="store_true",
        default=False,
        help="transparent background (default: False)",
    )
    parser.add_argument(
        "--title",
        required=False,
        type=str,
        default="panel",
        help="plot title (default: panel)",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        required=True,
        type=str,
        help="output file, panel.svg",
    )
    return parser
