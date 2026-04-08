"""Integration test for the full HATCHet3 Snakemake pipeline."""

import os
from pathlib import Path

import yaml
import pandas as pd
import pytest

from snakemake.api import SnakemakeApi
from snakemake.settings.types import (
    ConfigSettings,
    ResourceSettings,
)

from conftest import CBC_AVAILABLE
from simulate_bb_dir import simulate_bb_dir

pytestmark = pytest.mark.skipif(not CBC_AVAILABLE, reason="CBC solver not available")

SNAKEFILE = Path(__file__).resolve().parent.parent / "Snakefile"


@pytest.fixture(scope="module")
def snakemake_result(tmp_path_factory):
    """Run the full Snakemake pipeline on synthetic data."""
    base = tmp_path_factory.mktemp("snakemake")
    bb_dir = str(base / "bb_dir")
    genome_sizes = str(base / "genome.sizes")
    regions_bed = str(base / "regions.bed")
    workdir = base / "output"
    workdir.mkdir(exist_ok=True)

    ground_truth = simulate_bb_dir(bb_dir, genome_sizes, regions_bed)

    config = {
        "bb_dir": bb_dir,
        "genome_size": genome_sizes,
        "region_bed": regions_bed,
        "bbc_dir": "bbc",
        "result_dir": "results",
        "plot_dir": "summary",
        "log_dir": "logs",
        "verbosity": 1,
        "threads": 1,
        "cluster_bins": {
            "minK": 3,
            "maxK": 5,
            "t": 1e-6,
            "restarts": 3,
            "top_restarts": 2,
            "niters": 5,
            "min_tau": 50,
            "max_tau": 500,
            "baf_eps": 1e-3,
            "min_covar": 1e-3,
            "ig_alpha": 10.0,
            "tau_iters": 1,
            "log_rdr": False,
            "decode_method": "viterbi",
            "score_method": "icl",
            "init_method": "cna_plus_plus",
            "bal_lrt_alpha": 0.05,
            "bal_lrt_margin": 0.03,
            "filter_std": 2.0,
            "min_nbins": 10,
            "ub_nbins": 50,
        },
        "compute_cn": {
            "k": None,
            "mode": "ilp",
            "solver": "cbc",
            "timelimit": 60,
            "minClone": 2,
            "maxClone": 2,
            "diploid": True,
            "tetraploid": False,
            "segment": False,
            "reg_term": "MAXCN",
            "reg_steps": 3,
            "reg_bound": 0.15,
            "no_ampdel": False,
            "num_cnstates": -1,
            "diploidcmax": 6,
            "tetraploidcmax": 12,
            "min_prop": 0.01,
            "mrca": False,
            "max_degree": 3,
            "zero_cn_thres": 0.001,
            "purities": None,
            "fix_cn_dip": None,
            "fix_cn_tet": None,
            "cd_niters": 10,
            "cd_convergence_iters": 2,
            "cd_nseeds": 400,
            "cd_njobs": 1,
            "cd_seed": 42,
            "cd_tol": 0.001,
            "u_init": "dirichlet",
            "u_dir_alpha": 0.3,
            "solver_threads": None,
            "pool_size": 1,
            "pool_gap": None,
            "fcn_ci_alpha": 0.05,
            "min_ci_margin": 0.1,
        },
        "plot_cn": {
            "dpi": 100,
            "img_type": "png",
            "transparent": False,
            "keep_gap": False,
            "tail_alpha": 0.8,
            "center_alpha": 1.0,
            "onetail_area": 0.025,
            "maxlim_fcn": 30,
        },
    }

    config_file = base / "config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config, f)

    error = None
    try:
        with SnakemakeApi() as api:
            workflow = api.workflow(
                resource_settings=ResourceSettings(cores=1),
                config_settings=ConfigSettings(
                    configfiles=[config_file],
                ),
                snakefile=SNAKEFILE,
                workdir=workdir,
            )
            dag = workflow.dag()
            dag.execute_workflow()
    except Exception as exc:
        error = exc

    return str(workdir), ground_truth, error


class TestSnakemakePipeline:
    """Verify the full Snakemake pipeline runs and produces correct outputs."""

    def test_snakemake_exit_code(self, snakemake_result):
        workdir, _, error = snakemake_result
        assert error is None, f"Snakemake failed:\n{error}"

    def test_bbc_output_exists(self, snakemake_result):
        workdir, _, error = snakemake_result
        if error is not None:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "bbc", "bulk.bbc"))
        assert os.path.isfile(os.path.join(workdir, "bbc", "bulk.seg"))

    def test_best_ucn_exists(self, snakemake_result):
        workdir, _, error = snakemake_result
        if error is not None:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "results", "best.bbc.ucn"))
        assert os.path.isfile(os.path.join(workdir, "results", "best.seg.ucn"))

    def test_gammas_exists(self, snakemake_result):
        workdir, _, error = snakemake_result
        if error is not None:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "results", "gammas.tsv"))

    def test_ucn_format(self, snakemake_result):
        workdir, _, error = snakemake_result
        if error is not None:
            pytest.skip("Snakemake failed")
        bbc = pd.read_table(os.path.join(workdir, "results", "best.bbc.ucn"), sep="\t")
        assert "cn_normal" in bbc.columns
        assert "cn_clone1" in bbc.columns
        assert (bbc["cn_normal"] == "1|1").all()

    def test_logs_created(self, snakemake_result):
        workdir, _, error = snakemake_result
        if error is not None:
            pytest.skip("Snakemake failed")
        assert os.path.isfile(os.path.join(workdir, "logs", "cluster_bins.log"))
        assert os.path.isfile(os.path.join(workdir, "logs", "compute_cn.log"))
