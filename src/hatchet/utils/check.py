import sys
import os
import os.path
import subprocess
import textwrap
import importlib
from contextlib import contextmanager
from hatchet import config
from hatchet.utils.check_solver import main as check_solver


# from https://stackoverflow.com/questions/2125702
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        #sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            #sys.stderr = old_stderr


def _check_cmd(exe, *args, exe_with_path=None):
    # This function should never raise Exceptions unless it's a genuine implementation bug
    if exe_with_path is not None and exe_with_path.strip() != '':
        exe = exe_with_path
    cmd = [exe, *args]
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True)
        p.communicate()
        # We would like to verify the returncode, but certain programs like picard/tabix never return a 0!
        # assert p.returncode == 0
    except Exception:  # yes, a broad exception; See above
        return False
    else:
        return True


def _check_python_import(which):
    try:
        importlib.import_module(which)
    except ImportError:
        return False
    else:
        return True


_picard_jar_path = config.paths.picard or 'picard.jar'


# <HATCHet_command> => [(<dependency_name>, <error_message>, <boolean_func>, <func_args>..), ..] mapping
with suppress_stdout():
    CHECKS = {
        'count-reads': [
            (
                'tabix',
                'Please install tabix executable and either ensure its on your PATH, or its location specified in '
                'hatchet.ini as config.paths.tabix, or its location specified using the environment variable '
                'HATCHET_PATHS_TABIX',
                _check_cmd('tabix', exe_with_path=config.paths.tabix)
            ),
            (
                'mosdepth',
                'Please install mosdepth executable and either ensure its on your PATH, or its location specified in '
                'hatchet.ini as config.paths.mosdepth, or its location specified using the environment variable '
                'HATCHET_PATHS_MOSDEPTH',
                _check_cmd('mosdepth', '--version', exe_with_path=config.paths.mosdepth)
            )
        ],

        'phase-snps': [
            (
                'picard',
                'Please install picard.jar and ensure that its location is specified in hatchet.ini as '
                'config.paths.picard, or its location specified using the environment variable HATCHET_PATHS_PICARD. '
                'Also make sure "java" is on your path',
                os.path.exists(_picard_jar_path) and _check_cmd('java', '-jar', _picard_jar_path)
            ),
            (
                'shapeit',
                'Please install shapeit executable and either ensure its on your PATH, or its location specified in '
                'hatchet.ini as config.paths.shapeit, or its location specified using the environment variable '
                'HATCHET_PATHS_SHAPEIT',
                _check_cmd('shapeit', '--version', exe_with_path=config.paths.shapeit)
            )
        ],

        'compute-cn': [
            (
                'solver',
                'See http://compbio.cs.brown.edu/hatchet/README.html#using-a-solver',
                check_solver() is None
            )
        ],

        'cluster-bins': [
            (
                'hmmlearn',
                'Please install hmmlearn using pip/conda.',
                _check_python_import('hmmlearn')
            )
        ],

    }


def main(hatchet_cmds=None):
    all_ok = True
    hatchet_cmds = hatchet_cmds or CHECKS.keys()
    print('======================\nRunning HATCHet checks\n======================')
    for hatchet_cmd in hatchet_cmds:
        if hatchet_cmd in CHECKS:
            checks = CHECKS[hatchet_cmd]
            print(f'----------------------\nCommand: {hatchet_cmd}\n----------------------')
            for (cmd_name, msg, pred) in checks:
                if pred:
                    print(f'  {cmd_name} SUCCESSFUL')
                else:
                    msg = textwrap.fill(msg, initial_indent='    ', subsequent_indent='    ')
                    print(f'  {cmd_name} FAILED\n{msg}')
                    all_ok = False

    sys.exit(0 if all_ok else 1)


if __name__ == '__main__':
    main()
