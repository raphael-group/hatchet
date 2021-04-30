"""
The HATCHet package can be run as a module by invoking it as:
python -m hatchet <command> <arguments> ..

The commands currently supported are:

binBAM - description here
SNPCaller -
deBAF -
comBBo -
cluBB -
solve -
BBeval -
"""

import sys
import os.path
import hatchet
from hatchet.utils.binBAM import main as binBAM
from hatchet.utils.SNPCaller import main as SNPCaller
from hatchet.utils.deBAF import main as deBAF
from hatchet.utils.comBBo import main as comBBo
from hatchet.utils.cluBB import main as cluBB
from hatchet.utils.BBot import main as BBot
from hatchet.bin.HATCHet import main as solve
from hatchet.utils.BBeval import main as BBeval
from hatchet.utils.run import main as run

solve_bin = os.path.join(os.path.dirname(hatchet.__file__), 'solve')

commands = ('binBAM', 'SNPCaller', 'deBAF', 'comBBo', 'cluBB', 'BBot', 'solve', 'BBeval', 'run')


def print_usage():
    print('HATCHet v' + hatchet.__version__)
    print('Usage: hatchet <command> <arguments ..>')
    print('\nThe following commands are supported:\n ' + '\n '.join(commands))


def main():
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(0)

    command = sys.argv[1]
    args = sys.argv[2:]
    if command not in commands:
        print_usage()
        sys.exit(1)

    if command != 'solve':
        globals()[command](args)
    else:
        solve([solve_bin] + args)


if __name__ == '__main__':
    sys.exit(main())
