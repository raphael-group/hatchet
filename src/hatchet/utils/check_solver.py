import os.path
import tempfile
from importlib.resources import path
import traceback

from hatchet import config
import hatchet.data
from hatchet.bin.HATCHet import main as hatchet_main
from hatchet.utils.Supporting import log


this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def main(args=None):
    log(msg=f'# The solver you are currently using is {config.compute_cn.solver}\n')

    try:
        with tempfile.TemporaryDirectory() as tempdirname:
            with path(hatchet.data, 'sample.bbc') as bbc_path:
                input_files_prefix = os.path.splitext(bbc_path)[0]
                hatchet_main(
                    args=[
                        '-x',
                        os.path.join(tempdirname),
                        '-i',
                        input_files_prefix,
                        '-n4',
                        '-p',
                        '5',
                        '-v',
                        '3',
                        '-u',
                        '0',
                        '--mode',
                        '0',
                        '-r',
                        '6700',
                        '-j',
                        '1',
                        '-eD',
                        '6',
                        '-eT',
                        '12',
                        '-g',
                        '0.35',
                        '-l',
                        '0.6',
                        '-E',
                        '-bD',
                        '12',
                        '--uniqueclones',
                        '-P',
                        '0.6 0.2 0',
                    ]
                )

    except Exception as e:
        # write the exception message to a log file named check_solver.log in the current directory
        with open('check_solver.log', 'w') as f:
            f.write(traceback.print_stack())
            f.write('\n')
            f.write(str(e))
        return False
    else:
        log(msg=f'# Your current solver {config.compute_cn.solver} seems to be working correctly\n')
        return True


if __name__ == '__main__':
    main()
