# All supported HATCHet commands
commands = (
    'count-reads',
    'count-reads-fw',
    'genotype-snps',
    'count-alleles',
    'combine-counts',
    'combine-counts-fw',
    'cluster-bins',
    'plot-bins',
    'compute-cn',
    'plot-cn',
    'run',
    'check',
    'download-panel',
    'phase-snps',
    'cluster-bins-gmm',
    'plot-cn-1d2d',
    'plot-bins-1d2d',
)


# Support for old command names as they've been used in earlier versions of HATCHet
command_aliases = {
    'binBAM': 'count-reads-fw',
    'SNPCaller': 'genotype-snps',
    'deBAF': 'count-alleles',
    'comBBo': 'combine-counts-fw',
    'cluBB': 'cluster-bins-gmm',
    'BBot': 'plot-bins',
    'solve': 'compute-cn',
    'BBeval': 'plot-cn',
    'PhasePrep': 'download-panel',
    'Phase': 'phase-snps',
    'check-solver': 'check',
}
