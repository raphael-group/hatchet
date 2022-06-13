import sys
import os
import os.path
import subprocess
import textwrap
import importlib
from contextlib import contextmanager
import tempfile
import shutil

from hatchet import config
import hatchet.data
from hatchet.utils.check_solver import main as check_solver


# from https://stackoverflow.com/questions/2125702
@contextmanager
def suppress_stdout():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


def _check_cmd(exe_path, exe_name, *args):
    # This function should never raise Exceptions unless it's a genuine implementation bug
    # Only use exe and args that return a return code of 0
    # Use exe_path as '' if you expect <exe_name> to be on PATH
    cmd = [os.path.join(exe_path, exe_name), *args]
    try:
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        p.communicate()
        assert p.returncode == 0
    except Exception as e:
        return False
    else:
        return True


# Most command-line commands can be checked using _check_cmd(<folder_with_command>, <command>, '--version')
# Others, like below, need special handling because they have no simple invocations that return 0
def _check_tabix():
    with tempfile.TemporaryDirectory() as tempdirname:
        with importlib.resources.path(
            hatchet.data, 'sample.sorted.gff.gz'
        ) as gz_path:
            _temp_gz_path = os.path.join(tempdirname, 'sample.sorted.gff.gz')
            shutil.copy(gz_path, _temp_gz_path)
            return _check_cmd(
                config.paths.tabix, 'tabix', '-p', 'gff', _temp_gz_path, '-f'
            )


def _check_picard():
    picard_dir = config.paths.picard
    picard_java_flags = config.phase_snps.picard_java_flags
    picard_bin_path = os.path.join(picard_dir, 'picard')
    if os.path.exists(picard_bin_path):
        exe_path = f'{picard_dir}'
        exe_name = 'picard'
        args_pre = (picard_java_flags,)
    else:
        exe_path = ''  # Assume java is on PATH
        exe_name = 'java'
        args_pre = (
            picard_java_flags,
            '-jar',
            os.path.join(picard_dir, 'picard.jar'),
        )

    with tempfile.TemporaryDirectory() as tempdirname:
        with importlib.resources.path(
            hatchet.data, 'sample.sorted.bam'
        ) as bam_path:
            return _check_cmd(
                exe_path,
                exe_name,
                *args_pre,
                'BuildBamIndex',
                '--INPUT',
                bam_path,
                '--OUTPUT',
                f'{tempdirname}/sample.sorted.bam.bai',
            )


def _check_python_import(which):
    try:
        importlib.import_module(which)
    except ImportError:
        return False
    else:
        return True


# <HATCHet_command> => [(<dependency_name>, <success_message>, <failure_message>, <boolean_func>, <func_args>..), ..] mapping
CHECKS = {
    'count-reads': [
        (
            'tabix',
            '',
            'Please install tabix executable and either ensure its on your PATH, or its location specified in '
            'hatchet.ini as config.paths.tabix, or its location specified using the environment variable '
            'HATCHET_PATHS_TABIX',
            _check_tabix,
        ),
        (
            'mosdepth',
            '',
            'Please install mosdepth executable and either ensure its on your PATH, or its location specified in '
            'hatchet.ini as config.paths.mosdepth, or its location specified using the environment variable '
            'HATCHET_PATHS_MOSDEPTH',
            _check_cmd,
            config.paths.mosdepth,
            'mosdepth',
            '--version',
        ),
    ],
    'phase-snps': [
        (
            'picard',
            '',
            'Please install picard.jar and ensure that its location is specified in hatchet.ini as '
            'config.paths.picard, or its location specified using the environment variable HATCHET_PATHS_PICARD. '
            'Also make sure "java" is on your path',
            _check_picard,
        ),
        (
            'shapeit',
            '',
            'Please install shapeit executable and either ensure its on your PATH, or its location specified in '
            'hatchet.ini as config.paths.shapeit, or its location specified using the environment variable '
            'HATCHET_PATHS_SHAPEIT',
            _check_cmd,
            config.paths.shapeit,
            'shapeit',
            '--version',
        ),
    ],
    'compute-cn': [
        (
            'solver',
            f'Your selected solver "{config.compute_cn.solver}" seems to be working correctly',
            'See http://compbio.cs.brown.edu/hatchet/README.html#using-a-solver',
            check_solver,
        )
    ],
}


def main(hatchet_cmds=None):
    all_ok = True
    hatchet_cmds = hatchet_cmds or CHECKS.keys()
    print(
        '======================\nRunning HATCHet checks\n======================'
    )
    for hatchet_cmd in hatchet_cmds:
        if hatchet_cmd in CHECKS:
            checks = CHECKS[hatchet_cmd]
            print(
                f'----------------------\nCommand: {hatchet_cmd}\n----------------------'
            )
            for check in checks:
                cmd_name, success_msg, err_msg, func, args = (
                    check[0],
                    check[1],
                    check[2],
                    check[3],
                    check[4:],
                )
                with suppress_stdout():
                    pred = func(*args)
                if pred:
                    print(f'  {cmd_name} check SUCCESSFUL. {success_msg}')
                else:
                    msg = textwrap.fill(
                        err_msg,
                        initial_indent='    ',
                        subsequent_indent='    ',
                    )
                    print(f'  {cmd_name} check FAILED.\n{msg}')
                    all_ok = False

    sys.exit(0 if all_ok else 1)


if __name__ == '__main__':
    main()
