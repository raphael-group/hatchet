import os.path
import tempfile
from importlib.resources import path

from hatchet import config
import hatchet.data
from hatchet.bin.HATCHet import main as hatchet_main
from hatchet.utils.Supporting import log


this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')
solve_binary = os.path.join(os.path.dirname(hatchet.__file__), 'solve')


def main(args=None):
    log(msg=f'# The solver you are currently using is {config.compute_cn.solver}\n')

    with tempfile.TemporaryDirectory() as tempdirname:
        with path(hatchet.data, 'sample.bbc') as bbc_path:
            input_files_prefix = os.path.splitext(bbc_path)[0]
            hatchet_main(args=[
                solve_binary,
                '-x', os.path.join(tempdirname),
                '-i', input_files_prefix,
                '-n2',
                '-p', '5',
                '-v', '3',
                '-u', '0.03',
                '--mode', '1',
                '-r', '6700',
                '-j', '1',
                '-eD', '6',
                '-eT', '12',
                '-g', '0.35',
                '-l', '0.6'
            ])

    log(msg=f'# Your current solver {config.compute_cn.solver} seems to be working correctly\n')


if __name__ == '__main__':
    main()
