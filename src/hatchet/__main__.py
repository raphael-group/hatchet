"""
The HATCHet package can be run as a module by invoking it as:
python -m hatchet <command> <arguments> ..
"""

import sys
import os.path
import warnings
import hatchet

from hatchet.utils.preprocess import main as preprocess
from hatchet.utils.cluster_bins_kde import main as cluster_bins_kde
from hatchet.utils.count_reads import main as count_reads
from hatchet.utils.count_reads import main as count_reads_fw
from hatchet.utils.genotype_snps import main as genotype_snps
from hatchet.utils.count_alleles import main as count_alleles
from hatchet.utils.combine_counts import main as combine_counts
from hatchet.utils.combine_counts import main as combine_counts_fw
from hatchet.utils.cluster_bins import main as cluster_bins
from hatchet.utils.plot_bins import main as plot_bins
from hatchet.bin.HATCHet import main as compute_cn
from hatchet.utils.plot_cn import main as plot_cn
from hatchet.utils.run import main as run
from hatchet.utils.PhasePrep import main as PhasePrep
from hatchet.utils.Phase import main as Phase

solve_bin = os.path.join(os.path.dirname(hatchet.__file__), 'solve')

commands = ('count-reads', 'count-reads-fw', 'genotype-snps', 'count-alleles', 'combine-counts',
            'combine-counts-fw', 'cluster-bins', 'plot-bins','compute-cn', 'plot-cn', 'run',
            'PhasePrep', 'Phase', 'cluster-bins-kde')

def print_usage():
    print('HATCHet v' + hatchet.__version__)
    print('Usage: hatchet <command> <arguments ..>')
    print('\nThe following commands are supported:\n ' + '\n '.join(commands))


def main():
    # support for old command names as they've been used in earlier versions of HATCHet
    aliases = {
        'binBAM': 'count-reads-fw',
        'SNPCaller': 'genotype-snps',
        'deBAF': 'count-alleles',
        'comBBo': 'combine-counts-fw',
        'cluBB': 'cluster-bins',
        'BBot': 'plot-bins',
        'solve': 'compute-cn',
        'BBeval': 'plot-cn'
    }

    if len(sys.argv) < 2:
        print_usage()
        sys.exit(0)

    command = sys.argv[1]
    args = sys.argv[2:]

    if command in aliases:
        msg = f'The HATCHet command "{command}" has been replaced by "{aliases[command]}" and will be absent in ' \
              f'future releases. Please update your scripts accordingly.'
        warnings.warn(msg, FutureWarning)
    elif command not in commands:
        print_usage()
        sys.exit(1)

    command = aliases.get(command, command)
    if command != 'compute-cn':
        command = command.replace('-', '_')
        globals()[command](args)
    else:
        compute_cn([solve_bin] + args)


if __name__ == '__main__':
    sys.exit(main())
