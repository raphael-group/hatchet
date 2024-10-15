import os
import subprocess
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal

import hatchet
from hatchet.bin.HATCHet import main as main
from hatchet import config
from hatchet.utils.solve import solver_available

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, "data", "fw")
SOLVE = os.path.join(os.path.dirname(hatchet.__file__), "solve")


@pytest.fixture(scope="module")
def output_folder():
    out = os.path.join(this_dir, "out", "solver")
    os.makedirs(out, exist_ok=True)
    return out


@pytest.mark.skipif(
    os.getenv("GRB_LICENSE_FILE") is None, reason="No license for Gurobi found"
)
def test_solve_binary():
    cmd = [SOLVE, os.path.join(DATA_FOLDER, "bulk")]
    # -M 1 for ILP solve mode
    cmd.extend(
        [_x for _x in "-f -e 6 -p 400 -u 0.03 -r 6700 -M 1 -v 2 -c 5:1:1 -n 2".split()]
    )

    print(" ".join(cmd))

    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    p.communicate()
    assert p.returncode == 0


@pytest.mark.skipif(
    os.getenv("GRB_LICENSE_FILE") is None, reason="No license for Gurobi found"
)
def test_solve_command_cpp(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = "cpp"

    main(
        args=[
            "-x",
            os.path.join(output_folder),
            "-i",
            os.path.join(this_dir, "data/fw/bulk"),
            "-n2",
            "-p",
            "400",
            "-v",
            "3",
            "-u",
            "0.03",
            "--mode",
            "1",
            "-r",
            "6700",
            "-j",
            "8",
            "-eD",
            "6",
            "-eT",
            "12",
            "-g",
            "0.35",
            "-l",
            "0.6",
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, "best.bbc.ucn"), sep="\t")
    df2 = pd.read_csv(os.path.join(this_dir, "data", "fw", "best.bbc.ucn"), sep="\t")
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value


@pytest.mark.skipif(
    not solver_available("gurobipy"),
    reason="gurobipy solver not available for pyomo",
)
def test_solve_command_gurobipy(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = "gurobipy"

    main(
        args=[
            "-x",
            os.path.join(output_folder),
            "-i",
            os.path.join(this_dir, "data/fw/bulk"),
            "-n2",
            "-p",
            "400",
            "-v",
            "3",
            "-u",
            "0.03",
            "--mode",
            "1",
            "-r",
            "6700",
            "-j",
            "8",
            "-eD",
            "6",
            "-eT",
            "12",
            "-g",
            "0.35",
            "-l",
            "0.6",
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, "best.bbc.ucn"), sep="\t")
    df2 = pd.read_csv(os.path.join(this_dir, "data", "fw", "best.bbc.ucn"), sep="\t")
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value


@pytest.mark.skipif(
    not solver_available("gurobi"),
    reason="gurobi solver not available for pyomo",
)
def test_solve_command_gurobi(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = "gurobi"

    main(
        args=[
            "-x",
            os.path.join(output_folder),
            "-i",
            os.path.join(this_dir, "data/fw/bulk"),
            "-n2",
            "-p",
            "400",
            "-v",
            "3",
            "-u",
            "0.03",
            "--mode",
            "1",
            "-r",
            "6700",
            "-j",
            "8",
            "-eD",
            "6",
            "-eT",
            "12",
            "-g",
            "0.35",
            "-l",
            "0.6",
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, "best.bbc.ucn"), sep="\t")
    df2 = pd.read_csv(os.path.join(this_dir, "data", "fw", "best.bbc.ucn"), sep="\t")
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value


@pytest.mark.skipif(
    not solver_available("cbc"), reason="cbc solver not available for pyomo"
)
def test_solve_command_cbc(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = "cbc"

    main(
        args=[
            "-x",
            os.path.join(output_folder),
            "-i",
            os.path.join(this_dir, "data/fw/bulk"),
            "-n2",
            "-p",
            "400",
            "-v",
            "3",
            "-u",
            "0.03",
            "--mode",
            "1",
            "-r",
            "6700",
            "-j",
            "8",
            "-eD",
            "6",
            "-eT",
            "12",
            "-g",
            "0.35",
            "-l",
            "0.6",
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, "best.bbc.ucn"), sep="\t")
    df2 = pd.read_csv(os.path.join(this_dir, "data", "fw", "best.bbc.ucn"), sep="\t")
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value
