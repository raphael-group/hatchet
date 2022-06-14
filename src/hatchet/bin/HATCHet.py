import os
import sys
import argparse
import shutil
import subprocess
import shlex
from collections import Counter

import hatchet
from hatchet import config, __version__
from hatchet.utils.Supporting import ensure
from hatchet.utils.solve import solve
from hatchet.utils.solve.utils import segmentation


def parsing_arguments(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = ''
    parser = argparse.ArgumentParser(
        prog='hatchet compute-cn',
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        'SOLVER', type=str, nargs='?', help='Path to the executable solver'
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        help='Prefix path to seg and bbc input files (required)',
    )
    parser.add_argument(
        '-x',
        '--runningdir',
        type=str,
        required=False,
        default=config.compute_cn.runningdir,
        help='Running directory (default: ./)',
    )
    parser.add_argument(
        '-n',
        '--clones',
        type=str,
        required=False,
        default=config.compute_cn.clones,
        help='Either an estimated number of clones or an interval where the\nnumber of clones should be looked for given in the form LOWER,UPPER where LOWER and UPPER are two integer defining the interval (default: 2,8)',
    )
    parser.add_argument(
        '-f',
        '--noampdel',
        action='store_true',
        default=config.compute_cn.noampdel,
        required=False,
        help='Remove amp-del assumption where each mutated allele of every segment can be either amplified or deleted in all tumor clones w.r.t. base (2 for diploid and 4 for tetraploid) (default: use assumption)',
    )
    parser.add_argument(
        '-c',
        '--clonal',
        type=str,
        required=False,
        default=config.compute_cn.clonal,
        help='Clonal clusters to fix for tetraploid (default: automatically inferred)',
    )
    parser.add_argument(
        '-d',
        '--cnstates',
        type=int,
        required=False,
        default=config.compute_cn.cnstates,
        help='Maximum number of distinct copy-number states for each segment (default: None, no limit)',
    )
    parser.add_argument(
        '-eD',
        '--diploidcmax',
        type=int,
        required=False,
        default=config.compute_cn.diploidcmax,
        help='Maximum copy-number value overall segments (default: 6, 0 means inferred from scaled fractional copy numbers)',
    )
    parser.add_argument(
        '-eT',
        '--tetraploidcmax',
        type=int,
        required=False,
        default=config.compute_cn.tetraploidcmax,
        help='Maximum copy-number value overall segments (default: 12, 0 means inferred from scaled fractional copy numbers)',
    )
    parser.add_argument(
        '-ts',
        '--minsize',
        type=float,
        required=False,
        default=config.compute_cn.minsize,
        help='The minimum proportion of covered genome for potential clonal clusters (default: 0.008)',
    )
    parser.add_argument(
        '-tc',
        '--minchrs',
        type=int,
        required=False,
        default=config.compute_cn.minchrs,
        help='The minimum number of covered chromosomes for potential clonal clusters (default: 1)',
    )
    parser.add_argument(
        '-td',
        '--maxneutralshift',
        type=float,
        required=False,
        default=config.compute_cn.maxneutralshift,
        help='Maximum BAF shift for neutral cluster used to automatically infer the diploid/tetraploid cluster (default: 0.1)',
    )
    parser.add_argument(
        '--merge',
        action='store_true',
        default=config.compute_cn.merge,
        required=False,
        help='Merge the clusters (default: false)',
    )
    parser.add_argument(
        '-mR',
        '--mergeRDR',
        type=float,
        required=False,
        default=config.compute_cn.mergerdr,
        help='RDR tolerance used for finding the clonal copy numbers (default: 0.08)',
    )
    parser.add_argument(
        '-mB',
        '--mergeBAF',
        type=float,
        required=False,
        default=config.compute_cn.mergebaf,
        help='BAF tolerance used for finding the clonal copy numbers (default: 0.04)',
    )
    parser.add_argument(
        '-l',
        '--limitinc',
        type=float,
        required=False,
        default=config.compute_cn.limitinc,
        help='Upper bound to the relative increase of objective function. When there are significant small CNAs, their effect on the objective function may be confounded by only larger events, use this value to limit the relative increase of OBJ so that fitting small CNAs is more considered (default: None)',
    )
    parser.add_argument(
        '-g',
        '--ghostprop',
        type=float,
        required=False,
        default=config.compute_cn.ghostprop,
        help='Increasing proportion used to compute the value of the first ghost point added in the solution selection (default: 0.3)',
    )
    parser.add_argument(
        '-tR',
        '--toleranceRDR',
        type=float,
        required=False,
        default=config.compute_cn.tolerancerdr,
        help='RDR tolerance used for finding the clonal copy numbers (default: 0.08)',
    )
    parser.add_argument(
        '-tB',
        '--toleranceBAF',
        type=float,
        required=False,
        default=config.compute_cn.tolerancebaf,
        help='BAF tolerance used for finding the clonal copy numbers (default: 0.04)',
    )
    parser.add_argument(
        '-p',
        '--seeds',
        type=int,
        required=False,
        default=config.compute_cn.seeds,
        help='Number of seeds for coordinate-descent method (default: 400)',
    )
    parser.add_argument(
        '-j',
        '--jobs',
        type=int,
        required=False,
        default=config.compute_cn.jobs,
        help='Number of parallel jobs (default: maximum available on the machine)',
    )
    parser.add_argument(
        '-r',
        '--randomseed',
        type=int,
        required=False,
        default=config.compute_cn.randomseed,
        help='Random seed (default: None)',
    )
    parser.add_argument(
        '-s',
        '--timelimit',
        type=int,
        required=False,
        default=config.compute_cn.timelimit,
        help='Time limit for each ILP run (default: None)',
    )
    parser.add_argument(
        '-m',
        '--memlimit',
        type=int,
        required=False,
        default=config.compute_cn.memlimit,
        help='Memory limit for each ILP run (default: None)',
    )
    parser.add_argument(
        '-u',
        '--minprop',
        type=float,
        required=False,
        default=config.compute_cn.minprop,
        help='Minimum clone proporion in each sample (default: 0.03)',
    )
    parser.add_argument(
        '--maxiterations',
        type=int,
        required=False,
        default=config.compute_cn.maxiterations,
        help='Maximum number of iterations composed of C-step/U-step for each seed (default: 10)',
    )
    parser.add_argument(
        '--mode',
        type=int,
        required=False,
        default=config.compute_cn.mode,
        help='Solving mode among: Coordinate Descent + exact ILP (0), exact ILP only (1), and Coordinate-descent only (2) (default: 2)',
    )
    parser.add_argument(
        '--diploid',
        action='store_true',
        default=config.compute_cn.diploid,
        required=False,
        help='Force the tumor clones to be diploid without WGD (default: false)',
    )
    parser.add_argument(
        '--tetraploid',
        action='store_true',
        default=config.compute_cn.tetraploid,
        required=False,
        help='Force the tumor clones to be tetraploid with an occured WGD (default: false)',
    )
    parser.add_argument(
        '-v',
        '--verbosity',
        type=int,
        required=False,
        default=config.compute_cn.verbosity,
        help='Level of verbosity among: none (0), essential (1), verbose (2), and debug (3) (default: 1)',
    )
    parser.add_argument(
        '-V', '--version', action='version', version=f'%(prog)s {__version__}'
    )
    parser.add_argument(
        '-b',
        '--binwise',
        action='store_true',
        default=config.compute_cn.binwise,
        required=False,
        help="Use bin-wise objective function which requires more variables and constraints but accounts for cluster variances (default False). Only works with non-'cpp' solvers.",
    )
    args = parser.parse_args(args)

    if config.compute_cn.solver == 'cpp':
        if args.SOLVER is None:
            args.SOLVER = os.path.join(
                os.path.dirname(hatchet.__file__), 'solve'
            )
        if not os.path.isfile(args.SOLVER):
            raise ValueError(
                error('Solver not found in {}!'.format(args.SOLVER))
            )

    if config.compute_cn.solver == 'cpp' and args.binwise:
        raise ValueError(
            error(
                "The bin-wise objective is not supported for the solver 'cpp'. Please use a pyomo solver."
            )
        )

    clusters = args.input + '.seg'
    if not os.path.isfile(clusters):
        raise ValueError(error('SEG file {} not found!'.format(clusters)))
    bbc = args.input + '.bbc'
    if not os.path.isfile(bbc):
        raise ValueError(error('BBC file {} not found!'.format(bbc)))

    ln = -1
    un = -1
    if args.clones.isdigit():
        ln = int(args.clones)
        un = int(args.clones)
    else:
        parsed = args.clones.split(',')
        if (
            len(parsed) != 2
            or not parsed[0].isdigit()
            or not parsed[1].isdigit()
        ):
            raise ValueError(
                error('Wrong format for interval of clone numbers!')
            )
        else:
            ln = int(parsed[0])
            un = int(parsed[1])
    if ln < 2:
        raise ValueError(
            error(
                'The lower bound in the number of clones must be greater or equal than 2!'
            )
        )

    if args.cnstates is not None and args.cnstates <= 0:
        raise ValueError(
            error(
                'The maximum number of copy-number states should be default None or a positive non-zero integer!'
            )
        )
    if args.diploidcmax == 0:
        args.diploidcmax = None
    if args.diploidcmax is not None and args.diploidcmax <= 0:
        raise ValueError(
            error(
                'The maximum diploid copy number a positive non-zero integer!'
            )
        )
    if args.tetraploidcmax == 0:
        args.tetraploidcmax = None
    if args.tetraploidcmax is not None and args.tetraploidcmax <= 0:
        raise ValueError(
            error(
                'The maximum tetraploid copy number a positive non-zero integer!'
            )
        )
    if args.minsize < 0.0 or args.minsize > 1.0:
        raise ValueError(
            error(
                'The genome-size proportions for potential clonal clusters must be in [0, 1]!'
            )
        )
    if args.minchrs < 0 or args.minchrs > 22:
        raise ValueError(
            error(
                'The number of chromosomes for potential clonal clusters must be in {0, .., 22}!'
            )
        )
    if args.maxneutralshift < 0.0 or args.maxneutralshift > 1.0:
        raise ValueError(
            error(
                'The maximum BAF shift for neutral cluster must be in [0, 1]!'
            )
        )
    if args.toleranceRDR < 0.0 or args.toleranceRDR > 1.0:
        raise ValueError(
            error(
                'The RDR tolerance for finding clonal copy numbers must be in [0, 1]!'
            )
        )
    if args.toleranceBAF < 0.0 or args.toleranceBAF > 1.0:
        raise ValueError(
            error(
                'The BAF tolerance for finding clonal copy numbers must be in [0, 1]!'
            )
        )
    if args.mergeRDR < 0.0 or args.mergeRDR > 1.0:
        raise ValueError(
            error('The RDR tolerance for merging clusters must be in [0, 1]!')
        )
    if args.mergeBAF < 0.0 or args.mergeBAF > 1.0:
        raise ValueError(
            error('The BAF tolerance for merging clusters must be in [0, 1]!')
        )
    if args.limitinc is not None and (
        args.limitinc < 0.0 or args.limitinc > 1.0
    ):
        raise ValueError(error('The increasing limit must be in [0, 1]!'))
    if args.ghostprop < 0.0 or args.ghostprop > 1.0:
        raise ValueError(
            error(
                'The increasing proportion of the ghost point must be in [0, 1]!'
            )
        )
    if args.seeds <= 0:
        raise ValueError(
            error('The number of seeds should be a positive non-zero integer')
        )
    if args.jobs is not None and args.jobs <= 0:
        raise ValueError(
            error('The number of jobs should be a positive non-zero integer')
        )
    if args.randomseed is not None and args.randomseed <= 0:
        raise ValueError(
            error(
                'The random seed should be a positive non-zero integer or default None!'
            )
        )
    if args.timelimit is not None and args.timelimit <= 0:
        raise ValueError(
            error(
                'The time limit should be a positive non-zero integer or default None!'
            )
        )
    if args.memlimit is not None and args.memlimit <= 0:
        raise ValueError(
            error(
                'The memory limit should be a positive non-zero integer or default None!'
            )
        )
    if args.minprop is not None and (args.minprop < 0 or args.minprop > 0.3):
        raise ValueError(
            error(
                'The minimum proportion of clones on each sample must be in [0, 0.3]'
            )
        )
    if args.maxiterations is not None and args.maxiterations < 0:
        raise ValueError(
            error('The max-iteration number must be a positive integer!')
        )
    ensure(
        args.mode in (None, 0, 1, 2), 'The mode integer must be in (0, 1, 2)!'
    )

    if args.diploid and args.tetraploid:
        raise ValueError(
            error('Diploid and tetraploid cannot be forced at the same time!')
        )

    if not os.path.isdir(args.runningdir):
        raise ValueError(error('Running directory not found!'))
    if not 0 <= args.verbosity <= 3:
        raise ValueError(
            error('The verbosity level must be a value within 0,1,2,3!')
        )

    return {
        'solver': args.SOLVER,
        'input': args.input,
        'seg': clusters,
        'bbc': bbc,
        'ln': ln,
        'un': un,
        'clonal': args.clonal,
        'ampdel': not args.noampdel,
        'd': args.cnstates,
        'eD': args.diploidcmax,
        'eT': args.tetraploidcmax,
        'ts': args.minsize,
        'tc': args.minchrs,
        'td': args.maxneutralshift,
        'tR': args.toleranceRDR,
        'tB': args.toleranceBAF,
        'mR': args.mergeRDR,
        'mB': args.mergeBAF,
        'limit': args.limitinc,
        'g': args.ghostprop,
        'p': args.seeds,
        'j': args.jobs,
        'r': args.randomseed,
        's': args.timelimit,
        'm': args.memlimit,
        'u': args.minprop,
        'f': args.maxiterations,
        'M': args.mode,
        'x': args.runningdir,
        'diploid': args.diploid,
        'tetraploid': args.tetraploid,
        'v': args.verbosity,
        'binwise': args.binwise,
    }


def logArgs(args, width):
    text = '\n'
    for key in args:
        text += '\t{}: {}\n'.format(key, args[key])
    sys.stderr.write(log(text))


def main(args=None):
    sys.stderr.write(log('# Checking and parsing input arguments\n'))
    args = parsing_arguments(args)
    if args['v'] >= 2:
        sys.stdout.write(
            info(
                '\n'.join(
                    ['## {}:\t{}'.format(key, args[key]) for key in args]
                )
                + '\n'
            )
        )
    logArgs(args, 80)

    sys.stderr.write(log('# Reading and parsing clustered bins in BBC file\n'))
    bbc, bsamples = readBBC(args['bbc'])

    sys.stderr.write(log('# Reading and parsing bin clusters in SEG file\n'))
    seg, ssamples = readSEG(args['seg'])

    assert bsamples == ssamples, error(
        'Samples in BBC files does not match the ones in SEG file!'
    )
    samples = ssamples

    sys.stderr.write(log('# Computing the cluster sizes\n'))
    size = computeSizes(seg=seg, bbc=bbc, samples=samples)

    sys.stderr.write(
        log(
            '# Filtering clusters based on covered genome size and covered chromosomes\n'
        )
    )
    fbbc, fseg = filtering(
        bbc=bbc,
        seg=seg,
        size=size,
        ts=args['ts'],
        tc=args['tc'],
        mB=args['mB'],
        mR=args['mR'],
        samples=samples,
        v=args['v'],
    )

    sys.stderr.write(log('# Finding the neutral diploid/tetraploid cluster\n'))
    neutral = findNeutralCluster(
        seg=fseg, size=size, td=args['td'], samples=samples, v=args['v']
    )

    if not args['diploid'] and not args['tetraploid']:
        sys.stderr.write(log('# Running diploid\n'))
        diploidObjs = runningDiploid(neutral=neutral, args=args)

        if args['clonal'] is None:
            sys.stderr.write(
                log('# Finding clonal clusters and their copy numbers\n')
            )
            clonal, scale = findClonalClusters(
                fseg=fseg,
                neutral=neutral,
                size=size,
                tB=args['tB'],
                tR=args['tR'],
                samples=samples,
                v=args['v'],
            )
        else:
            sys.stderr.write(log('# Parsing given clonal copy numbers\n'))
            clonal, scale = parseClonalClusters(
                clonal=args['clonal'],
                fseg=fseg,
                neutral=neutral,
                size=size,
                samples=samples,
                v=args['v'],
            )

        if len(clonal) > 0:
            sys.stderr.write(log('# Running tetraploid\n'))
            tetraploidObjs = runningTetraploid(
                clonal=clonal, scale=scale, size=size, args=args
            )

            sys.stderr.write(log('# Selecting best solution\n'))
            select(
                diploid=diploidObjs,
                tetraploid=tetraploidObjs,
                v=args['v'],
                rundir=args['x'],
                g=args['g'],
                limit=args['limit'],
            )
        else:
            sys.stderr.write(
                warning(
                    '# No potential clonal patterns found, the input is likely to be diploid.\n If the heuristic failed to identify a clonal cluster due to high noisy in the data, there are two main parameters which user may tune:\n\t1. Increase the values of tB and tR which define the thresholds to predict the RDR and BAF of clonal clusters.\n\t2. Decrease the value of -ts which define the minimum coverage of the genome for clusters to be considered potentially clonal. As such, more clusters will be considered. \nLast, please assess the quality of the clustering through cluster-bins and increase the corresponding thresholds (-tR and -tB) to avoid overfitting of the data and overclustering.\n'
                )
            )

            sys.stderr.write(log('# Selecting best diploid solution\n'))
            selectDiploid(
                diploid=diploidObjs,
                v=args['v'],
                rundir=args['x'],
                g=args['g'],
                limit=args['limit'],
            )

    elif args['diploid']:
        sys.stderr.write(log('# Running diploid\n'))
        diploidObjs = runningDiploid(neutral=neutral, args=args)

        sys.stderr.write(log('# Selecting best diploid solution\n'))
        selectDiploid(
            diploid=diploidObjs,
            v=args['v'],
            rundir=args['x'],
            g=args['g'],
            limit=args['limit'],
        )

    elif args['tetraploid']:
        if args['clonal'] is None:
            sys.stderr.write(
                log('# Finding clonal clusters and their copy numbers\n')
            )
            clonal, scale = findClonalClusters(
                fseg=fseg,
                neutral=neutral,
                size=size,
                tB=args['tB'],
                tR=args['tR'],
                samples=samples,
                v=args['v'],
            )
        else:
            sys.stderr.write(log('# Parsing given clonal copy numbers\n'))
            clonal, scale = parseClonalClusters(
                clonal=args['clonal'],
                fseg=fseg,
                neutral=neutral,
                size=size,
                samples=samples,
                v=args['v'],
            )

        if len(clonal) > 0:
            sys.stderr.write(log('# Running tetraploid\n'))
            tetraploidObjs = runningTetraploid(
                clonal=clonal, scale=scale, size=size, args=args
            )

            sys.stderr.write(log('# Selecting best tetraploid solution\n'))
            selectTetraploid(
                tetraploid=tetraploidObjs,
                v=args['v'],
                rundir=args['x'],
                g=args['g'],
                limit=args['limit'],
            )
        else:
            sys.stderr.write(
                warning(
                    '# No potential clonal patterns found, the input is likely to be diploid.\n If the heuristic failed to identify a clonal cluster due to high noisy in the data, there are two main parameters which user may tune:\n\t1. Increase the values of tB and tR which define the thresholds to predict the RDR and BAF of clonal clusters.\n\t2. Decrease the value of -ts which define the minimum coverage of the genome for clusters to be considered potentially clonal. As such, more clusters will be considered. \nLast, please assess the quality of the clustering through cluster-bins and increase the corresponding thresholds (-tR and -tB) to avoid overfitting of the data and overclustering.\n'
                )
            )


def readBBC(filename):
    samples = set()
    bbc = {}
    with open(filename, 'r') as f:
        for line in f:
            parsed = line.strip().split()
            if len(parsed) == 11 and line[0] != '#':
                chro = str(parsed[0])
                start = int(parsed[1])
                end = int(parsed[2])
                sample = str(parsed[3])
                if chro not in bbc:
                    bbc[chro] = {}
                if (start, end) not in bbc[chro]:
                    bbc[chro][start, end] = {}
                if sample not in bbc[chro][start, end]:
                    bbc[chro][start, end][sample] = {}
                bbc[chro][start, end][sample]['rdr'] = float(parsed[4])
                bbc[chro][start, end][sample]['snps'] = int(parsed[5])
                bbc[chro][start, end][sample]['cov'] = float(parsed[6])
                bbc[chro][start, end][sample]['alpha'] = int(parsed[7])
                bbc[chro][start, end][sample]['beta'] = int(parsed[8])
                bbc[chro][start, end][sample]['baf'] = float(parsed[9])
                bbc[chro][start, end][sample]['cluster'] = str(parsed[10])
                samples.add(sample)
            elif line[0] != '#' and len(parsed) > 0:
                raise ValueError(
                    error('Found a bad-formatted line: "{}"'.format(line))
                )

    for c in bbc:
        for seg in bbc[c]:
            for p in samples:
                assert p in bbc[c][seg], error(
                    'ERROR: bin ({}:{}-{}) not defined in sample {}'.format(
                        c, seg[0], seg[1], p
                    )
                )

    return bbc, samples


def readSEG(filename):
    samples = set()
    seg = {}
    with open(filename, 'r') as f:
        for line in f:
            parsed = line.strip().split()
            if len(parsed) == 9 and line[0] != '#':
                idx = str(parsed[0])
                sample = str(parsed[1])
                if idx not in seg:
                    seg[idx] = {}
                if sample not in seg[idx]:
                    seg[idx][sample] = {}
                seg[idx][sample]['bins'] = int(parsed[2])
                seg[idx][sample]['rdr'] = float(parsed[3])
                seg[idx][sample]['snps'] = int(parsed[4])
                seg[idx][sample]['cov'] = float(parsed[5])
                seg[idx][sample]['alpha'] = int(parsed[6])
                seg[idx][sample]['beta'] = int(parsed[7])
                seg[idx][sample]['baf'] = float(parsed[8])
                samples.add(sample)
            elif line[0] != '#' and len(parsed) > 0:
                raise ValueError(
                    error('Found a bad-formatted line: "{}"'.format(line))
                )

    for i in seg:
        for p in samples:
            assert p in seg[i], error(
                'ERROR: cluster {} not defined in sample {}'.format(i, p)
            )

    return seg, samples


def computeSizes(seg, bbc, samples):
    sample = list(samples)[0]
    size = {
        idx: sum(
            float(b[1] - b[0])
            for c in bbc
            for b in bbc[c]
            if bbc[c][b][sample]['cluster'] == idx
        )
        for idx in seg
    }
    for idx in seg:
        for p in samples:
            seg[idx][p]['avgbaf'] = (
                sum(
                    bbc[c][b][p]['baf'] * float(b[1] - b[0])
                    for c in bbc
                    for b in bbc[c]
                    if bbc[c][b][p]['cluster'] == idx
                )
                / size[idx]
            )
    return size


def filtering(bbc, seg, size, ts, tc, mB, mR, samples, v):
    sample = list(samples)[0]
    totsize = float(sum(seg[1] - seg[0] for c in bbc for seg in bbc[c]))
    chrs = {
        idx: set(
            c
            for c in bbc
            for seg in bbc[c]
            if bbc[c][seg][sample]['cluster'] == idx
        )
        for idx in seg
    }
    assert sum(size[idx] for idx in size) == totsize
    printrdrbaf = lambda s: '-'.join(
        [
            '({}, {})'.format(seg[s][p]['rdr'], seg[s][p]['baf'])
            for p in samples
        ]
    )

    if v >= 3:
        sys.stderr.write(debug('### Clusters and their features:\n'))
        # sys.stderr.write(debug('\n'.join(['### {}: SIZE= {}\t#CHRS= {}\t(RDR, BAF)= {}'.format(idx, size[idx], len(chrs[idx]), printrdrbaf(idx)) for idx in seg]) + '\n'))

    selected = set(idx for idx in seg)
    merge_selected = set()
    for idx in sorted(selected, key=(lambda i: size[i]), reverse=True):
        merge = None
        for ins in sorted(
            merge_selected, key=(lambda i: size[i]), reverse=True
        ):
            if set(
                abs(seg[idx][p]['rdr'] - seg[ins][p]['rdr']) <= mR
                and abs(seg[idx][p]['baf'] - seg[ins][p]['baf']) <= mB
                for p in samples
            ) == {True}:
                merge = ins
                break
        if merge is None:
            merge_selected.add(idx)
            if v >= 3:
                sys.stderr.write(
                    debug(
                        '### {}: SIZE= {}\t#CHRS= {}\t(RDR, BAF)= {}\n'.format(
                            idx, size[idx], len(chrs[idx]), printrdrbaf(idx)
                        )
                    )
                )
        else:
            size[merge] += size[idx]
            chrs[merge] = chrs[merge] | chrs[idx]
            if v >= 3:
                sys.stderr.write(
                    warning(
                        '### {} merged with {}: SIZE= {}\t#CHRS= {}\t(RDR, BAF)= {}\n'.format(
                            idx,
                            merge,
                            size[idx],
                            len(chrs[idx]),
                            printrdrbaf(idx),
                        )
                    )
                )
    selected = merge_selected

    selected = set(
        idx
        for idx in selected
        if size[idx] >= totsize * ts and len(chrs[idx]) >= tc
    )

    fbbc = {
        c: {
            seg: bbc[c][seg]
            for seg in bbc[c]
            if bbc[c][seg][sample]['cluster'] in selected
        }
        for c in bbc
    }
    fseg = {idx: seg[idx] for idx in selected}

    if v >= 1:
        sys.stderr.write(
            info(
                '# Selected clusters: '
                + ', '.join([idx for idx in selected])
                + '\n'
            )
        )
    if v >= 2:
        sys.stderr.write(info('## Features of selected clusters:\n'))
        sys.stderr.write(
            info(
                '\n'.join(
                    [
                        '## {}: SIZE= {}\t#CHRS= {}\t(RDR, BAF)= {}'.format(
                            idx, size[idx], len(chrs[idx]), printrdrbaf(idx)
                        )
                        for idx in selected
                    ]
                )
                + '\n'
            )
        )

    return fbbc, fseg


def findNeutralCluster(seg, size, td, samples, v):
    selected = set(
        idx
        for idx in seg
        if set(p for p in samples if abs(0.5 - seg[idx][p]['baf']) <= td)
        == samples
    )
    if len(selected) == 0:
        raise ValueError(
            error(
                'None potential neutral cluster found with given parameters!'
            )
        )
    neutral = max(selected, key=(lambda k: size[k]))

    if v >= 2:
        sys.stderr.write(
            info(
                '## Cluster selected as neutral (diploid/tetraploid) is {}\n'.format(
                    neutral
                )
            )
        )

    return neutral


def runningDiploid(neutral, args):
    results = []
    args['c'] = '{}:1:1'.format(neutral)
    args['e'] = args['eD']
    basecmd = makeBaseCMD(args, args['eD']) + ' -c ' + args['c']

    for n in range(args['ln'], args['un'] + 1):
        if args['v'] >= 2:
            sys.stderr.write(
                log('## Running diploid with {} clones\n'.format(n))
            )
        outprefix = os.path.join(args['x'], 'results.diploid.n{}'.format(n))
        results.append((n, execute(args, basecmd, n, outprefix), outprefix))

    return results


def makeBaseCMD(args, e):
    basecmd = '{} {}'.format(args['solver'], args['input'])
    if args['ampdel']:
        basecmd += ' -f '
    if args['d'] is not None:
        basecmd += ' -d {}'.format(args['d'])
    if e is not None:
        basecmd += ' -e {}'.format(e)
    if args['j'] is not None:
        basecmd += ' -j {}'.format(args['j'])
    if args['p'] is not None:
        basecmd += ' -p {}'.format(args['p'])
    if args['u'] is not None:
        basecmd += ' -u {}'.format(args['u'])
    if args['m'] is not None:
        basecmd += ' -m {}'.format(args['m'])
    if args['s'] is not None:
        basecmd += ' -s {}'.format(args['s'])
    if args['f'] is not None:
        basecmd += ' -i {}'.format(args['f'])
    if args['r'] is not None:
        basecmd += ' -r {}'.format(args['r'])
    if args['M'] is not None:
        basecmd += ' -M {}'.format(args['M'])
    basecmd += ' -v 2'
    # if args['td'] is not None:
    #    basecmd += ' -t {}'.format(args['td'])
    return basecmd


def findClonalClusters(fseg, neutral, size, tB, tR, samples, v):
    topbaf = {
        p: float(min(fseg[idx][p]['avgbaf'] for idx in fseg)) for p in samples
    }
    heigh = {
        p: (abs(fseg[neutral][p]['avgbaf'] - topbaf[p]) / 4.0) > tB
        for p in samples
    }

    location = {}
    level = {}
    for idx in fseg:
        locp = []
        levp = {}

        for p in samples:
            t = topbaf[p]
            n = fseg[neutral][p]['avgbaf']
            x = fseg[idx][p]['avgbaf']
            d = {
                'top': abs(t - x),
                'mid': abs(((t + n) / 2.0) - x),
                'midbot': abs(((t + 3 * n) / 4.0) - x),
                'bot': abs(n - x),
            }

            if d['top'] < d['bot'] and abs(d['top'] - d['bot']) > tB:
                levp[p] = ['top']
            elif d['bot'] < d['top'] and abs(d['top'] - d['bot']) > tB:
                levp[p] = ['bot']
            else:
                levp[p] = ['top', 'bot']

            c = argmin(d)
            locp.append(c)
            if c != 'top' and d['top'] <= tB:
                locp.append('top')
            if c != 'mid' and d['mid'] <= tB:
                locp.append('mid')
            if c != 'midbot' and d['midbot'] <= tB:
                locp.append('midbot')
            if c != 'bot' and d['bot'] <= tB:
                locp.append('bot')

        count = Counter([levv for p in samples for levv in levp[p]])
        count = sorted(count.keys(), key=(lambda x: count[x]), reverse=True)

        for lev in count:
            if False not in set(lev in levp[p] for p in samples):
                level[idx] = lev
                locp = Counter(locp)
                loc = argmax(locp)
                location[idx] = loc
                break

    clusters = sorted(
        [idx for idx in fseg if idx != neutral and idx in location],
        key=(lambda i: size[i]),
        reverse=True,
    )
    allclonal = [
        (2, 0),
        (2, 1),
        (3, 2),
        (4, 2),
        (1, 0),
        (3, 0),
        (3, 1),
        (4, 0),
        (4, 1),
        (5, 0),
    ]
    found_pattern = []
    best_pattern = {}
    best_scale = ()
    best_value = 0

    for cluster in clusters:
        rightpos = sum(
            (fseg[cluster][p]['rdr'] - fseg[neutral][p]['rdr']) > tR
            for p in samples
        )
        leftpos = sum(
            (fseg[cluster][p]['rdr'] - fseg[neutral][p]['rdr']) < -tR
            for p in samples
        )

        eqbaf = (
            lambda p: abs(fseg[cluster][p]['baf'] - fseg[neutral][p]['baf'])
            <= tB
        )
        eqrdr = (
            lambda p: abs(fseg[cluster][p]['rdr'] - fseg[neutral][p]['rdr'])
            <= tR
        )

        if True in set(eqbaf(p) and heigh[p] for p in samples):
            continue
        if True in set(eqbaf(p) and eqrdr(p) and heigh[p] for p in samples):
            continue

        if rightpos == len(samples):
            if location[cluster] == 'bot' or location[cluster] == 'midbot':
                options = [(3, 2), (4, 2)]
            elif location[cluster] == 'mid' or location[cluster] == 'top':
                options = [(4, 2), (3, 2)]
            else:
                assert False
            # if sum(abs(fseg[cluster][p]['baf'] - fseg[neutral][p]['baf']) <= tB for p in samples) == len(samples):
            #     options = [(4, 4)] + options
        elif leftpos == len(samples):
            if level[cluster] == 'top' and location[cluster] == 'top':
                options = [(2, 0)]
            elif level[cluster] == 'top':
                options = [(2, 0), (2, 1)]
            elif level[cluster] == 'bot':
                options = [(2, 1)]
            else:
                assert False
            # if sum(abs(fseg[cluster][p]['baf'] - fseg[neutral][p]['baf']) <= tB for p in samples) == len(samples):
            #     options = [(1, 1), (0, 0)] + options
        else:
            options = []
            continue

        if v >= 3:
            sys.stderr.write(
                debug(
                    '### Potential clonal cluster {} with copy numbers {}\n'.format(
                        cluster, options
                    )
                )
            )

        regPurity = (
            lambda v: 1.0
            if 1.0 <= v <= 1.05
            else (0.0 if -0.05 <= v <= 0.0 else v)
        )
        calcPurity = lambda d, c, r: regPurity(
            float(2 * d - 2 * r) / float(2 * r + 2 * d - c * d)
        )
        calcScalingFactor = lambda p, d: float(2.0 + 2.0 * p) / float(d)
        calcFraction = lambda p, cn: float(2.0 * (1.0 - p) + sum(cn) * p)
        calcRDR = lambda p, cn, s: calcFraction(p, cn) / float(s)
        calcBAF = lambda p, cn: float(
            1.0 * (1.0 - p) + min(cn) * p
        ) / calcFraction(p, cn)

        for opt in options:
            purity = {
                p: calcPurity(
                    fseg[neutral][p]['rdr'], sum(opt), fseg[cluster][p]['rdr']
                )
                for p in samples
            }
            if False in set(0.0 <= purity[p] <= 1.0 for p in samples):
                continue
            scaling = {
                p: calcScalingFactor(purity[p], fseg[neutral][p]['rdr'])
                for p in samples
            }
            if False in set(scaling[p] >= 0.0 for p in samples):
                continue
            curr_pattern = {}
            curr_pattern[neutral] = (2, 2)
            curr_pattern[cluster] = opt
            curr_scale = (neutral, cluster)
            curr_value = size[neutral] + size[cluster]
            for cn in [a for a in allclonal if a != opt]:
                estRDR = {
                    p: calcRDR(purity[p], cn, scaling[p]) for p in samples
                }
                estBAF = {p: calcBAF(purity[p], cn) for p in samples}
                checkRDR = (
                    lambda r, i: set(
                        p
                        for p in samples
                        if abs(r[p] - fseg[i][p]['rdr']) <= tR
                    )
                    == samples
                )
                checkBAF = (
                    lambda b, i: set(
                        p
                        for p in samples
                        if abs(b[p] - fseg[i][p]['baf']) <= tB
                    )
                    == samples
                )
                candidates = [
                    idx
                    for idx in clusters
                    if idx not in curr_pattern
                    and checkRDR(estRDR, idx)
                    and checkBAF(estBAF, idx)
                ]
                if len(candidates) > 0:
                    choice = max(candidates, key=(lambda i: size[i]))
                    curr_pattern[choice] = cn
                    curr_value += size[choice]

            if v >= 2:
                if curr_pattern not in found_pattern:
                    sys.stderr.write(
                        info(
                            '## Found pattern of size {}: {}\n'.format(
                                curr_value, curr_pattern
                            )
                        )
                    )
            if curr_pattern not in found_pattern:
                found_pattern.append(curr_pattern)

            if curr_value > best_value:
                best_pattern = curr_pattern
                best_value = curr_value
                best_scale = curr_scale

    if v >= 1:
        sys.stderr.write(
            info(
                '## Chosen pattern of size {}: {}\n'.format(
                    best_value, best_pattern
                )
            )
        )

    return best_pattern, best_scale


def parseClonalClusters(clonal, fseg, neutral, size, samples, v):
    given = clonal.split(',')
    if len(given) < 2:
        raise RuntimeError(
            error('At least two clonal clusters must be provided!')
        )

    for e in given:
        p = e.split(':')
        if len(p) != 3:
            raise RuntimeError(error('Wrong format of clonal clusters!'))
        if p[0] not in set(fseg):
            raise RuntimeError(
                error(
                    'A specified clonal cluster does not exist or is not selected! {}'.format(
                        p[0]
                    )
                )
            )

    tmp = given[0].split(':')
    neutral = tmp[0]
    """
    if int(tmp[1]) != 2 or int(tmp[1]) != 2:
        raise RuntimeError(error('The first clonal cluster must be the neutral with copy numbers 2:2!'))
    """
    tmp = given[1].split(':')
    second = tmp[0]
    cn = int(tmp[1]) + int(tmp[2])

    calcPurity = lambda d, c, r: float(2 * d - 2 * r) / float(
        2 * r + 2 * d - c * d
    )
    calcScalingFactor = lambda p, d: float(2.0 + 2.0 * p) / float(d)

    purity = {
        p: calcPurity(fseg[neutral][p]['rdr'], cn, fseg[second][p]['rdr'])
        for p in samples
    }
    if False in set(0.0 <= purity[p] <= 1.0 for p in samples):
        raise RuntimeError(
            error(
                'The specified clonal clusters do not allow for scaling because resulting purity is {}!'.format(
                    purity
                )
            )
        )
    scaling = {
        p: calcScalingFactor(purity[p], fseg[neutral][p]['rdr'])
        for p in samples
    }
    if False in set(scaling[p] >= 0.0 for p in samples):
        raise RuntimeError(
            error(
                'The specified clonal clusters do not allow for scaling because resulting scaling factor is {}!'.format(
                    scaling
                )
            )
        )

    pattern = {
        e.split(':')[0]: (int(e.split(':')[1]), int(e.split(':')[2]))
        for e in given
    }
    scale = (neutral, second)

    return pattern, scale


def runningTetraploid(clonal, scale, size, args):
    results = []
    cn = '{}:{}:{}'.format(scale[0], clonal[scale[0]][0], clonal[scale[0]][1])
    cn += ',{}:{}:{}'.format(
        scale[1], clonal[scale[1]][0], clonal[scale[1]][1]
    )
    if len(clonal) > 2:
        cn = (
            cn
            + ','
            + ','.join(
                [
                    '{}:{}:{}'.format(s, clonal[s][0], clonal[s][1])
                    for s in clonal
                    if s not in scale
                ]
            )
        )

    args['c'] = cn
    args['e'] = args['eT']
    basecmd = makeBaseCMD(args, args['eT']) + ' -c {}'.format(cn)

    for n in range(args['ln'], args['un'] + 1):
        if args['v'] >= 2:
            sys.stderr.write(
                log('## Running tetraploid with {} clones\n'.format(n))
            )
        outprefix = os.path.join(args['x'], 'results.tetraploid.n{}'.format(n))
        results.append((n, execute(args, basecmd, n, outprefix), outprefix))

    return results


def execute_python(solver, args, n, outprefix):

    bbc_out_file = outprefix + '.bbc.ucn.tsv'
    seg_out_file = outprefix + '.seg.ucn.tsv'

    _mode = ('both', 'ilp', 'cd')[args['M']]

    obj, cA, cB, u, cluster_ids, sample_ids = solve(
        solver=solver,
        clonal=args['c'],
        bbc_file=args['bbc'],
        seg_file=args['seg'],
        n=n,
        solve_mode=_mode,
        d=-1 if args['d'] is None else args['d'],
        cn_max=args['e'],
        mu=args['u'],
        diploid_threshold=0.1,
        ampdel=args['ampdel'],
        n_seed=args['p'],
        n_worker=args['j'],
        random_seed=args['r'],
        max_iters=args['f'],
        timelimit=args['s'],
        binwise=args['binwise'],
    )

    segmentation(
        cA,
        cB,
        u,
        cluster_ids,
        sample_ids,
        bbc_file=args['bbc'],
        bbc_out_file=bbc_out_file,
        seg_out_file=seg_out_file,
    )

    if args['v'] >= 1:
        sys.stderr.write(
            info('# Best objective found with {} clones: {}\n'.format(n, obj))
        )

    return obj


def execute(args, basecmd, n, outprefix):
    if config.compute_cn.solver != 'cpp':
        return execute_python(config.compute_cn.solver, args, n, outprefix)
    progressbar = ProgressBar(
        total=args['p'], length=min(50, args['p']), verbose=False
    )
    cmd = basecmd + ' -n {} -o {}'.format(n, outprefix)
    if args['v'] >= 3:
        sys.stderr.write(debug('### Running command: {}\n'.format(cmd)))

    FNULL = open(os.devnull, 'w')
    process = subprocess.Popen(
        shlex.split(cmd),
        stdout=FNULL,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    buffer = []

    if args['v'] >= 2:
        while True:
            nextline = process.stderr.readline().strip()
            buffer.append(nextline)
            if ';' in nextline:
                progressbar.progress(
                    advance=True,
                    msg='Obj value {}'.format(
                        float(nextline.split()[-2][:-1])
                    ),
                )
            if not nextline and process.poll() is not None:
                break

    stdout, stderr = process.communicate()
    exitcode = process.returncode
    if exitcode == 0:
        obj = -1
        for line in reversed(buffer):
            if 'Final objective' in line:
                for e in line.split():
                    if isfloat(e):
                        obj = float(e)
                break
        if obj >= 0:
            if args['v'] >= 1:
                sys.stderr.write(
                    info(
                        '# Best objective found with {} clones: {}\n'.format(
                            n, obj
                        )
                    )
                )
        else:
            raise RuntimeError(
                error(
                    'Failed to parse the output of the following command because the final objective was not found: \n\t\t{}\n'.format(
                        cmd
                    )
                )
            )
    else:
        if any('GRBException' in line for line in buffer):
            msg = '\nYou likely have a licensing issue with Gurobi. Please run `hatchet check-solver` to ensure that the solver is working correctly.'
        else:
            msg = '\nUnexpected error during solve. Please run `hatchet check-solver` to ensure that the solver is working correctly.'
        raise RuntimeError(
            error(
                'The following command failed: \n\t\t{}\nwith {}\n{}'.format(
                    cmd, buffer, msg
                )
            )
        )

    return obj


def select(diploid, tetraploid, v, rundir, g, limit):
    assert len(diploid) == len(tetraploid), error(
        'The number of diploid and tetraploid results must be the same'
    )
    dscores = {}
    tscores = {}

    if len(diploid) == 1 or len(diploid) == 2:
        for dip in diploid:
            dscores[dip[0]] = dip[1]
        if v >= 2:
            sys.stderr.write(
                info(
                    '## Objective value is used as scores for diploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Diploid with {} clones - OBJ: {} - score: {}'.format(
                                d[0], d[1], dscores[d[0]]
                            )
                            for d in diploid
                        ]
                    )
                    + '\n'
                )
            )

        for tet in tetraploid:
            tscores[tet[0]] = tet[1]
        if v >= 2:
            sys.stderr.write(
                info(
                    '## Objective value is used as scores for tetraploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Tetraploid with {} clones - OBJ: {} - score: {}'.format(
                                t[0], t[1], tscores[t[0]]
                            )
                            for t in tetraploid
                        ]
                    )
                    + '\n'
                )
            )
    else:
        for i, dip in enumerate(diploid):
            if i == 0:
                dscores[dip[0]] = forward(diploid, i, g, limit)
            elif i == len(diploid) - 1:
                dscores[dip[0]] = backward(diploid, i, g, limit)
            else:
                dscores[dip[0]] = central(diploid, i, g, limit)

        if v >= 2:
            sys.stderr.write(
                info(
                    '## Scores approximating second derivative for diploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Diploid with {} clones - OBJ: {} - score: {}'.format(
                                d[0], d[1], dscores[d[0]]
                            )
                            for d in diploid
                        ]
                    )
                    + '\n'
                )
            )

        for i, tet in enumerate(tetraploid):
            if i == 0:
                tscores[tet[0]] = estimate_forward(
                    tetraploid, i, g, diploid, limit
                )
            elif i == len(tetraploid) - 1:
                tscores[tet[0]] = backward(tetraploid, i, g, limit)
            else:
                tscores[tet[0]] = central(tetraploid, i, g, limit)

        if v >= 2:
            sys.stderr.write(
                info(
                    '## Scores approximating second derivative for tetraploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Tetraploid with {} clones - OBJ: {} - score: {}'.format(
                                t[0], t[1], tscores[t[0]]
                            )
                            for t in tetraploid
                        ]
                    )
                    + '\n'
                )
            )

    dchosen = max(diploid, key=(lambda x: dscores[x[0]]))
    dbout = os.path.join(rundir, 'chosen.diploid.bbc.ucn')
    shutil.copy2(dchosen[2] + '.bbc.ucn.tsv', dbout)
    dsout = os.path.join(rundir, 'chosen.diploid.seg.ucn')
    shutil.copy2(dchosen[2] + '.seg.ucn.tsv', dsout)
    if v >= 1:
        sys.stderr.write(
            info(
                '# The chosen diploid solution has {} clones with OBJ: {} and score: {}\n'.format(
                    dchosen[0], dchosen[1], dscores[dchosen[0]]
                )
            )
        )
    if v >= 2:
        sys.stderr.write(
            info(
                '## The related-diploid resulting files are copied to {} and {}\n'.format(
                    dbout, dsout
                )
            )
        )

    tchosen = max(tetraploid, key=(lambda x: tscores[x[0]]))
    tbout = os.path.join(rundir, 'chosen.tetraploid.bbc.ucn')
    shutil.copy2(tchosen[2] + '.bbc.ucn.tsv', tbout)
    tsout = os.path.join(rundir, 'chosen.tetraploid.seg.ucn')
    shutil.copy2(tchosen[2] + '.seg.ucn.tsv', tsout)
    if v >= 1:
        sys.stderr.write(
            info(
                '# The chosen tetraploid solution has {} clones with OBJ: {} and score: {}\n'.format(
                    tchosen[0], tchosen[1], tscores[tchosen[0]]
                )
            )
        )
    if v >= 2:
        sys.stderr.write(
            info(
                '## The related-tetraploid resulting files are copied to {} and {}\n'.format(
                    tbout, tsout
                )
            )
        )

    bbest = os.path.join(rundir, 'best.bbc.ucn')
    sbest = os.path.join(rundir, 'best.seg.ucn')
    if tchosen[0] < dchosen[0]:
        shutil.copy2(tchosen[2] + '.bbc.ucn.tsv', bbest)
        shutil.copy2(tchosen[2] + '.seg.ucn.tsv', sbest)
        sys.stderr.write(
            info(
                '# The chosen solution is tetraploid with {} clones and is written in {} and {}\n'.format(
                    tchosen[0], bbest, sbest
                )
            )
        )
    #        if tchosen[0] == tetraploid[0][0] and False not in set(tscores[tet[0]] <= 0.05 for i, tet in enumerate(tetraploid) if i > 0):
    #            sys.stderr.write(warning("### HATCHet inferred the presence of a single tetraploid tumor clone and there are two possible cases for this inference:\n\t1. This is indeed the true clonal composition of the tumor. This case can be typically confirmed by the fact that \n\t\t\t(I) The chosen number of clones for diploid solutions is also pretty small (3-4 clones)\n\t\t\t(II) the objective value for tetraploid solutions keeps to be slightly decreasing when increasing the number of clones \n\t\t\t(III) the objective value of tetraploid solutions is not hugely different from the corresponding objective value of diploid solutions with the same number of clones. \n\t2. The heuristic to identify a clonal cluster failed. When this is the case, there are three typical observations: \n\t\t\t(I) The chosen number of clones for diploid solutions can be higher\n\t\t\t(II) the objective value does not almost vary (or only very slightly) when increasing the number of clones \n\t\t\t(III) the objective value is typically much higher than the corresponding values for diploid solutions.\n\nWhen the solution with a single tetraploid tumor clone appears suspicious, please use the CBB plot generated by plot-bins to evaluate the clonal cluster idetified by the heuristic. If the identified clonal cluster is clearly wrong there are two options to fix the inference:\n\t1. Use the parameters which control the heuristic of this method, more specifically increase/decrease the values of the thresholds tR and tB to consider higher/lower noise in the data, or decrease/increase the minimum coverage of the genome -ts used to select potential clonal clusters in order to have less/more clonal clusters to consider. Typically, adding more cluster helps the heuristic to have more information available.\n\t2. Improve the clustering of cluster-bins, e.g. increase the values of the thresholds tR and tB of cluster-bins to avoid overfitting and overclustering."))
    else:
        shutil.copy2(dchosen[2] + '.bbc.ucn.tsv', bbest)
        shutil.copy2(dchosen[2] + '.seg.ucn.tsv', sbest)
        sys.stderr.write(
            info(
                '# The chosen solution is diploid with {} clones and is written in {} and {}\n'.format(
                    dchosen[0], bbest, sbest
                )
            )
        )


def selectDiploid(diploid, v, rundir, g, limit):
    dscores = {}

    if len(diploid) == 1 or len(diploid) == 2:
        for dip in diploid:
            dscores[dip[0]] = dip[1]
        if v >= 2:
            sys.stderr.write(
                info(
                    '## Objective value is used as scores for diploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Diploid with {} clones - OBJ: {} - score: {}'.format(
                                d[0], d[1], dscores[d[0]]
                            )
                            for d in diploid
                        ]
                    )
                    + '\n'
                )
            )
    else:
        for i, dip in enumerate(diploid):
            if i == 0:
                dscores[dip[0]] = forward(diploid, i, g, limit)
            elif i == len(diploid) - 1:
                dscores[dip[0]] = backward(diploid, i, g, limit)
            else:
                dscores[dip[0]] = central(diploid, i, g, limit)

        if v >= 2:
            sys.stderr.write(
                info(
                    '## Scores approximating second derivative for diploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Diploid with {} clones - OBJ: {} - score: {}'.format(
                                d[0], d[1], dscores[d[0]]
                            )
                            for d in diploid
                        ]
                    )
                    + '\n'
                )
            )

    dchosen = max(diploid, key=(lambda x: dscores[x[0]]))
    dbout = os.path.join(rundir, 'chosen.diploid.bbc.ucn')
    shutil.copy2(dchosen[2] + '.bbc.ucn.tsv', dbout)
    dsout = os.path.join(rundir, 'chosen.diploid.seg.ucn')
    shutil.copy2(dchosen[2] + '.seg.ucn.tsv', dsout)
    if v >= 1:
        sys.stderr.write(
            info(
                '# The chosen diploid solution has {} clones with OBJ: {} and score: {}\n'.format(
                    dchosen[0], dchosen[1], dscores[dchosen[0]]
                )
            )
        )
    if v >= 2:
        sys.stderr.write(
            info(
                '## The related-diploid resulting files are copied to {} and {}\n'.format(
                    dbout, dsout
                )
            )
        )

    bbest = os.path.join(rundir, 'best.bbc.ucn')
    sbest = os.path.join(rundir, 'best.seg.ucn')

    shutil.copy2(dchosen[2] + '.bbc.ucn.tsv', bbest)
    shutil.copy2(dchosen[2] + '.seg.ucn.tsv', sbest)
    sys.stderr.write(
        info(
            '# The chosen solution is diploid with {} clones and is written in {} and {}\n'.format(
                dchosen[0], bbest, sbest
            )
        )
    )


def selectTetraploid(tetraploid, v, rundir, g, limit):
    tscores = {}

    if len(tetraploid) == 1 or len(tetraploid) == 2:
        for tet in tetraploid:
            tscores[tet[0]] = tet[1]
        if v >= 2:
            sys.stderr.write(
                info(
                    '## Objective value is used as scores for tetraploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Tetraploid with {} clones - OBJ: {} - score: {}'.format(
                                t[0], t[1], tscores[t[0]]
                            )
                            for t in tetraploid
                        ]
                    )
                    + '\n'
                )
            )
    else:
        for i, tet in enumerate(tetraploid):
            if i == 0:
                tscores[tet[0]] = forward(tetraploid, i, g, limit)
            elif i == len(tetraploid) - 1:
                tscores[tet[0]] = backward(tetraploid, i, g, limit)
            else:
                tscores[tet[0]] = central(tetraploid, i, g, limit)

        if v >= 2:
            sys.stderr.write(
                info(
                    '## Scores approximating second derivative for tetraploid results\n'
                )
            )
            sys.stderr.write(
                info(
                    '\n'.join(
                        [
                            '## Tetraploid with {} clones - OBJ: {} - score: {}'.format(
                                t[0], t[1], tscores[t[0]]
                            )
                            for t in tetraploid
                        ]
                    )
                    + '\n'
                )
            )

    tchosen = max(tetraploid, key=(lambda x: tscores[x[0]]))
    tbout = os.path.join(rundir, 'chosen.tetraploid.bbc.ucn')
    shutil.copy2(tchosen[2] + '.bbc.ucn.tsv', tbout)
    tsout = os.path.join(rundir, 'chosen.tetraploid.seg.ucn')
    shutil.copy2(tchosen[2] + '.seg.ucn.tsv', tsout)
    if v >= 1:
        sys.stderr.write(
            info(
                '# The chosen tetraploid solution has {} clones with OBJ: {} and score: {}\n'.format(
                    tchosen[0], tchosen[1], tscores[tchosen[0]]
                )
            )
        )
    if v >= 2:
        sys.stderr.write(
            info(
                '## The related-tetraploid resulting files are copied to {} and {}\n'.format(
                    tbout, tsout
                )
            )
        )

    bbest = os.path.join(rundir, 'best.bbc.ucn')
    sbest = os.path.join(rundir, 'best.seg.ucn')
    shutil.copy2(tchosen[2] + '.bbc.ucn.tsv', bbest)
    shutil.copy2(tchosen[2] + '.seg.ucn.tsv', sbest)
    sys.stderr.write(
        info(
            '# The chosen solution is tetraploid with {} clones and is written in {} and {}\n'.format(
                tchosen[0], bbest, sbest
            )
        )
    )


def safediv(v):
    return v if v > 0 else 1.0


def forward(f, i, g, limit):
    if limit is None:
        left = g
    else:
        left = min(g, limit)
    right = float(max(f[i][1] - f[i + 1][1], 0.0) / safediv(f[i][1]))
    return left - right


def estimate_forward(f, i, g, knw, limit):
    left = max(
        g, float(max(knw[i][1] - knw[i + 1][1], 0.0) / safediv(knw[i][1]))
    )
    if limit is not None:
        left = min(left, limit)
    right = float(max(f[i][1] - f[i + 1][1], 0.0) / safediv(f[i][1]))
    return left - right


def central(f, i, g, limit):
    if limit is None:
        left = float(max(f[i - 1][1] - f[i][1], 0.0) / safediv(f[i - 1][1]))
    else:
        left = min(
            limit,
            float(max(f[i - 1][1] - f[i][1], 0.0) / safediv(f[i - 1][1])),
        )
    right = float(max(f[i][1] - f[i + 1][1], 0.0) / safediv(f[i][1]))
    return left - right


def backward(f, i, g, limit):
    if limit is None:
        left = float(max(f[i - 1][1] - f[i][1], 0.0) / safediv(f[i - 1][1]))
    else:
        left = min(
            limit,
            float(max(f[i - 1][1] - f[i][1], 0.0) / safediv(f[i - 1][1])),
        )
    right = g
    return left - right


def argmin(d):
    return min(d, key=d.get)


def argmax(d):
    return max(d, key=d.get)


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def error(msg):
    return '{}{}{}'.format('\033[91m\033[1m', msg, '\033[0m')


def log(msg):
    return '{}{}{}'.format('\033[95m\033[1m', msg, '\033[0m')


def warning(msg):
    return '{}{}{}'.format('\033[93m', msg, '\033[0m')


def info(msg):
    return '{}{}{}'.format('\033[96m', msg, '\033[0m')


def debug(msg):
    return '{}{}{}'.format('\033[92m', msg, '\033[0m')


class ProgressBar:
    def __init__(
        self,
        total,
        length,
        counter=0,
        verbose=False,
        decimals=1,
        fill=chr(9608),
        prefix='Progress:',
        suffix='Complete',
    ):
        self.total = total
        self.length = length
        self.decimals = decimals
        self.fill = fill
        self.prefix = prefix
        self.suffix = suffix
        self.counter = counter
        self.verbose = verbose

    def progress(self, advance=True, msg=''):
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            self.counter += 1
        percent = ('{0:.' + str(self.decimals) + 'f}').format(
            100 * (self.counter / float(self.total))
        )
        filledLength = int(self.length * self.counter // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = '%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix)
        if not self.verbose:
            toprint = rewind + result + ' [%s]' % (msg)
        else:
            toprint = rewind + msg + '\n' + result
        write(toprint)
        flush()
        if self.counter == self.total:
            write('\n')
            flush()


if __name__ == '__main__':
    main()
