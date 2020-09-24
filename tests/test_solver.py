import os
import subprocess
import pytest
import hatchet

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')
SOLVE = os.path.join(os.path.dirname(hatchet.__file__), 'solve')


@pytest.mark.skipif(os.getenv('GRB_LICENSE_FILE') is None, reason='No license for Gurobi found')
def test_solver():
    cmd = [
        SOLVE,
        os.path.join(DATA_FOLDER, 'bulk')
    ]
    # -M 1 for ILP solve mode
    cmd.extend([_x for _x in '-f -e 6 -p 400 -u 0.03 -r 6700 -M 1 -v 2 -c 1:1:1 -n 2'.split()])

    print ' '.join(cmd)

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()
    assert p.returncode == 0
