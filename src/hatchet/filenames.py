"""HATCHet input/output constant filenames and subdirectories.

Producers and consumers of a file must import the SAME name from here, so a
rename cannot silently break a cross-stage contract. Fixed names are module
constants; names that embed run parameters (ploidy, clone count, K, sample) are
helper functions.

Notes:
    Purely user-named outputs (e.g. plot out_prefix) are not centralized here.
    The legacy ``evaluate_pool_solutions`` reader references filenames with no
    current producer and is intentionally not represented here.
"""

# =============================================================================
# Subdirectories
# =============================================================================
LABELS_DIR = "labels"  # cluster-bins per-K BBC/SEG
PLOTS_DIR = "plots"  # cluster-bins + compute-cn plots
TRACES_DIR = "traces"  # cluster-bins EM traces
INIT_DIAG_DIR = "init_diag"  # cluster-bins init diagnostics
SOLS_DIR = "sols"  # compute-cn per-(ploidy, n) solutions


# =============================================================================
# cluster-bins inputs (bb_dir)
# =============================================================================
BB_TSV_GZ = "bb.tsv.gz"
SAMPLE_IDS = "sample_ids.tsv"
BB_RDR_NPZ = "bb.rdr.npz"
BB_DEPTH_NPZ = "bb.depth.npz"
BB_A_ALLELE_NPZ = "bb.Aallele.npz"
BB_B_ALLELE_NPZ = "bb.Ballele.npz"
BB_T_ALLELE_NPZ = "bb.Tallele.npz"


# =============================================================================
# cluster-bins outputs (bbc_dir)
# =============================================================================
BULK_BBC = "bulk.bbc"
BULK_SEG = "bulk.seg"
BB_PHASED_TSV_GZ = "bb.phased.tsv.gz"
MODEL_SCORES_TSV = "model_scores.tsv"
MODEL_SCORES_PDF = "model_scores.pdf"
ELBO_TRACES_PDF = "elbo_traces.pdf"
HMM_INIT_PDF = "hmm_init.pdf"


def bulk_k_bbc(k) -> str:
    """Per-K clustered BBC under labels/."""
    return f"bulk{k}.bbc"


def bulk_k_seg(k) -> str:
    """Per-K segment SEG under labels/."""
    return f"bulk{k}.seg"


def bulk_k_phased(k) -> str:
    """Per-K phased bins under labels/."""
    return f"bulk{k}.bb.phased.tsv.gz"


def k_em_trace(k) -> str:
    """Per-K EM parameter trace under traces/."""
    return f"K{k}.em_trace.npz"


def k_plot(k) -> str:
    """Per-K 1D/2D plot under plots/."""
    return f"K{k}.pdf"


def bulk_k_plot(k) -> str:
    """Best-K plot copied to the top-level bbc_dir."""
    return f"bulk.K{k}.pdf"


def init_pdf(name) -> str:
    """Per-init-method HMM init diagnostic under init_diag/."""
    return f"{name}_init.pdf"


# =============================================================================
# compute-cn outputs (result_dir)
# =============================================================================
GAMMAS = "gammas.tsv"
SCALING_2D_PDF = "scaling_2d.pdf"
SUMMARY_TSV = "summary.tsv"
BEST_BBC_UCN = "best.bbc.ucn"
BEST_SEG_UCN = "best.seg.ucn"
MODEL_SELECTION_PDF = "model_selection.pdf"
POOL_PDF = "pool.pdf"
# sols/
OBJECTIVES_TSV = "objectives.tsv"
U0_SEEDS_TSV = "u0_seeds.tsv"


def solver_input(ploidy) -> str:
    """Per-ploidy solver input under sols/."""
    return f"solver_input.{ploidy}.tsv"


def results_bbc_ucn(ploidy, n) -> str:
    """Per-(ploidy, n) bin-level UCN."""
    return f"results.{ploidy}.n{n}.bbc.ucn.tsv"


def results_seg_ucn(ploidy, n) -> str:
    """Per-(ploidy, n) segment-level UCN."""
    return f"results.{ploidy}.n{n}.seg.ucn.tsv"


def chosen_bbc_ucn(ploidy) -> str:
    """Model-selected bin-level UCN for one ploidy."""
    return f"chosen.{ploidy}.bbc.ucn"


def chosen_seg_ucn(ploidy) -> str:
    """Model-selected segment-level UCN for one ploidy."""
    return f"chosen.{ploidy}.seg.ucn"


def ploidy_n_subdir(ploidy, n) -> str:
    """Per-(ploidy, n) leaf directory name (under sols/ and plots/)."""
    return f"{ploidy}_n{n}"


def solution_stem(solve_mode, sol_id) -> str:
    """Per-solution filename stem (extension appended by the caller: .tsv/.nwk/.json)."""
    return f"{solve_mode}_{sol_id}"


def solution_tsv(solve_mode, sol_id) -> str:
    """Per-solution detail TSV under a sols/<ploidy>_n<n>/ directory."""
    return solution_stem(solve_mode, sol_id) + ".tsv"


def pool_pdf(pid, ploidy, n) -> str:
    """Pool CNP panel for one (ploidy, n)."""
    return f"{pid}.pool_{ploidy}_n{n}.pdf"


# =============================================================================
# evaluate outputs (out_dir)
# =============================================================================
SOMATIC_SNVS_TSV = "somatic_snvs.tsv"
POOL_EVAL_TSV = "pool_eval.tsv"
EVAL_SUMMARY_TSV = "eval_summary.tsv"


def vaf_1d_pdf(sample) -> str:
    """Per-sample VAF 1D plot."""
    return f"{sample}.vaf_1d.pdf"


# =============================================================================
# logs
# =============================================================================
RUNTIME_LOG = "runtime.log"


def command_log(command) -> str:
    """Per-command file log written into the command's output directory."""
    return f"{command}.log"
