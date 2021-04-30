"""
The HATCHet package can be run as a module by invoking it as:
python -m hatchet <command> <arguments> ..
"""

import sys
import os.path
import hatchet
from hatchet.utils.count_reads import main as count_reads
from hatchet.utils.genotype_snps import main as genotype_snps
from hatchet.utils.count_alleles import main as count_alleles
from hatchet.utils.combine_counts import main as combine_counts
from hatchet.utils.cluster_bins import main as cluster_bins
from hatchet.utils.plot_bins import main as plot_bins
from hatchet.bin.HATCHet import main as compute_cn
from hatchet.utils.plot_cn import main as plot_cn

solve_bin = os.path.join(os.path.dirname(hatchet.__file__), 'solve')


if __name__ == '__main__':

    commands = ('count-reads', 'genotype-snps', 'count-alleles', 'combine-counts', 'cluster-bins', 'plot-bins',
                'compute-cn', 'plot-cn')
    if len(sys.argv) < 2:
        print('Usage: python -m hatchet <command> <arguments ..>')
        sys.exit(0)

    command = sys.argv[1]
    args = sys.argv[2:]
    if command not in commands:
        print('The following commands are supported: ' + ' '.join(commands))
        sys.exit(1)

    if command != 'compute_cn':
        command = command.replace('-', '_')
        globals()[command](args)
    else:
        compute_cn([solve_bin] + args)
