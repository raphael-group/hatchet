"""
The HATCHet package can be run as a module by invoking it as:
python -m hatchet <command> <arguments> ..
"""

import sys
import warnings
import hatchet

from hatchet.utils.commands import commands, command_aliases

from hatchet.utils.count_reads import main as count_reads  # noqa: F401
from hatchet.utils.count_reads_fw import main as count_reads_fw  # noqa: F401
from hatchet.utils.genotype_snps import main as genotype_snps  # noqa: F401
from hatchet.utils.count_alleles import main as count_alleles  # noqa: F401
from hatchet.utils.combine_counts import main as combine_counts  # noqa: F401
from hatchet.utils.combine_counts_fw import main as combine_counts_fw  # noqa: F401
from hatchet.utils.cluster_bins_gmm import main as cluster_bins_gmm  # noqa: F401
from hatchet.utils.cluster_bins import main as cluster_bins  # noqa: F401

from hatchet.utils.plot_bins import main as plot_bins  # noqa: F401
from hatchet.utils.plot_bins_1d2d import main as plot_bins_1d2d  # noqa: F401
from hatchet.bin.HATCHet import main as compute_cn  # noqa: F401
from hatchet.utils.plot_cn import main as plot_cn  # noqa: F401
from hatchet.utils.plot_cn_1d2d import main as plot_cn_1d2d  # noqa: F401

from hatchet.utils.run import main as run  # noqa: F401
from hatchet.utils.download_panel import main as download_panel  # noqa: F401
from hatchet.utils.phase_snps import main as phase_snps  # noqa: F401

from hatchet.utils.check import main as check  # noqa: F401


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

    if command in command_aliases:
        msg = (
            f'The HATCHet command "{command}" has been replaced by "{command_aliases[command]}" and will be absent in '
            f'future releases. Please update your scripts accordingly.'
        )
        warnings.warn(msg, FutureWarning)
    elif command not in commands:
        print_usage()
        sys.exit(1)

    command = command_aliases.get(command, command)
    command = command.replace('-', '_')
    globals()[command](args)


if __name__ == '__main__':
    sys.exit(main())
