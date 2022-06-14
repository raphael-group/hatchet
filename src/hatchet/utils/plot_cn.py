import sys
import os
import argparse
import shutil
import subprocess
import shlex
import math
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.colors as col
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter
from collections import deque
import itertools
from itertools import cycle

from hatchet.utils.Supporting import to_tuple
from hatchet import config, __version__

mpl.use('Agg')
plt.style.use('ggplot')
sns.set_style('whitegrid')


def parsing_arguments(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = ''
    parser = argparse.ArgumentParser(
        prog='hatchet plot-cn',
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        'INPUT', help='One or more space-separated files in CN_BBC format'
    )
    parser.add_argument(
        '-n',
        '--patientnames',
        required=False,
        default=config.plot_cn.patientnames,
        type=str,
        help='One or more space-separated patient names (default: inferred from filenames)',
    )
    parser.add_argument(
        '-u',
        '--minu',
        required=False,
        default=config.plot_cn.minu,
        type=float,
        help='Minimum proportion of a CNA to be considered subclonal (default: 0.2)"',
    )
    parser.add_argument(
        '-x',
        '--rundir',
        required=False,
        default=config.plot_cn.rundir,
        type=str,
        help='Running directory (default: current directory)',
    )
    parser.add_argument(
        '-b',
        '--baseCN',
        required=False,
        default=config.plot_cn.basecn,
        type=int,
        help='Base copy number (default: inferred from tumor ploidy)',
    )
    parser.add_argument(
        '-sC',
        '--figsizeclones',
        required=False,
        default=config.plot_cn.figsizeclones,
        type=str,
        help='Size of clone plots in the form "(X-SIZE, Y-SIZE)"',
    )
    parser.add_argument(
        '-sP',
        '--figsizecn',
        required=False,
        default=config.plot_cn.figsizecn,
        type=str,
        help='Size of CN plots in the form "(X-SIZE, Y-SIZE)"',
    )
    parser.add_argument(
        '-sG',
        '--figsizegrid',
        required=False,
        default=config.plot_cn.figsizegrid,
        type=str,
        help='Size of grid plots in the form "(X-SIZE, Y-SIZE)"',
    )
    parser.add_argument(
        '-rC',
        '--resolutionclones',
        required=False,
        default=config.plot_cn.resolutionclones,
        type=int,
        help='Number of bins to merge together for plotting clone profiles (default: 100)"',
    )
    parser.add_argument(
        '-rP',
        '--resolutioncn',
        required=False,
        default=config.plot_cn.resolutioncn,
        type=int,
        help='Number of bins to merge together for plotting proportions (default: 500)"',
    )
    parser.add_argument(
        '-rG',
        '--resolutiongrid',
        required=False,
        default=config.plot_cn.resolutiongrid,
        type=int,
        help='Number of bins to merge together in grids (default: 100)"',
    )
    parser.add_argument(
        '-e',
        '--threshold',
        required=False,
        default=config.plot_cn.threshold,
        type=float,
        help='Threshold used to classify a tumor into either diploid or tetraploid (default: 3.0)"',
    )
    parser.add_argument(
        '--ymax',
        required=False,
        default=config.plot_cn.ymax,
        type=int,
        help='Maximum values in y-axis (default: automatically inferred)"',
    )
    parser.add_argument(
        '--ymin',
        required=False,
        default=config.plot_cn.ymin,
        type=int,
        help='Minimum values in y-axis (default: automatically inferred)"',
    )
    parser.add_argument(
        '--clonepalette',
        required=False,
        default=config.plot_cn.clonepalette,
        type=str,
        help='Palette for coloring the clones among Set1, Set2, Set3, Paired (default: Set1)"',
    )
    parser.add_argument(
        '--linkage',
        required=False,
        default=config.plot_cn.linkage,
        type=str,
        help='Linkage method used for clustering (default: single, available (single, complete, average, weighted, centroid, median, ward) from SciPy)"',
    )
    parser.add_argument(
        '-V', '--version', action='version', version=f'%(prog)s {__version__}'
    )
    args = parser.parse_args(args)

    if len(args.INPUT.split()) == 0:
        raise ValueError(error('Please specify at least one sample as input!'))
    if args.patientnames is None:
        patientnames = {
            fil: os.path.basename(fil) for fil in args.INPUT.split()
        }
    else:
        patientnames = {
            f: n for f, n in zip(args.INPUT.split(), args.patientnames.split())
        }
    if len(args.INPUT.split()) != len(set(patientnames.values())):
        raise ValueError(
            error(
                'Multiple patients have the same name but they should unique!'
            )
        )
    if args.figsizeclones is not None:
        figsizeclones = to_tuple(
            args.figsizeclones, 'Wrong format of figsizeclones!'
        )
    if args.figsizecn is not None:
        figsizecn = to_tuple(args.figsizecn, 'Wrong format of figsizecn!')
    if args.figsizegrid is not None:
        figsizegrid = to_tuple(
            args.figsizegrid, 'Wrong format of figsizegrid!'
        )

    if not os.path.isdir(args.rundir):
        raise ValueError(error('Running directory does not exist!'))
    if not 0.0 <= args.minu <= 1.0:
        raise ValueError(
            error(
                'The minimum proportion for subclonal CNAs must be in [0, 1]!'
            )
        )
    if args.baseCN is not None and args.baseCN < 2:
        raise ValueError(error('Base CN must be greater or equal than 2!'))
    if args.resolutionclones is not None and args.resolutionclones < 1:
        raise ValueError(error('Resolution must be greater than 1!'))
    if args.resolutioncn is not None and args.resolutioncn < 1:
        raise ValueError(error('Resolution must be greater than 1!'))
    if args.resolutiongrid is not None and args.resolutiongrid < 1:
        raise ValueError(error('Resolution must be greater than 1!'))
    if args.threshold < 0:
        raise ValueError(error('Threshold must be positive!'))
    if args.linkage not in {
        'single',
        'complete',
        'average',
        'weighted',
        'centroid',
        'median',
        'ward',
    }:
        raise ValueError(error('Unknown linkage method!'))

    if args.clonepalette == 'Set1':
        pal = plt.cm.Set1
    elif args.clonepalette == 'Set2':
        pal = plt.cm.Set2
    elif args.clonepalette == 'Set3':
        pal = plt.cm.Set3
    elif args.clonepalette == 'Paired':
        pal = plt.cm.Paired
    else:
        raise ValueError(error('Unknown clone palette!'))

    return {
        'input': args.INPUT.split(),
        'names': patientnames,
        'rundir': args.rundir,
        'minu': args.minu,
        'base': args.baseCN,
        'clonefigsize': figsizeclones,
        'propsfigsize': figsizecn,
        'clusterfigsize': figsizegrid,
        'profileres': args.resolutionclones,
        'cnres': args.resolutioncn,
        'clusterres': args.resolutiongrid,
        'threshold': args.threshold,
        'linkage': args.linkage,
        'ymax': args.ymax,
        'ymin': args.ymin,
        'clonepalette': pal,
    }


def main(args=None):
    sys.stderr.write(log('# Checking and parsing input arguments\n'))
    args = parsing_arguments(args)
    sys.stdout.write(
        info(
            '\n'.join(['## {}:\t{}'.format(key, args[key]) for key in args])
            + '\n'
        )
    )

    sys.stderr.write(log('# Read BBC.UCN files\n'))
    tumors, clones, props = readUCN(args['input'], args['names'])

    sys.stderr.write(
        log(
            '# Compute purity and tumor ploidy of each sample from each patient\n'
        )
    )
    infbase = pp(tumors, clones, props, args)

    if args['base'] is None:
        sys.stderr.write(
            log('# The estimated basic copy number for each patient is\n')
        )
        sys.stderr.write(
            info(
                '\n'.join(['## {}: {}'.format(b, infbase[b]) for b in infbase])
                + '\n'
            )
        )
        base = infbase
    else:
        base = {pat: args['base'] for pat in tumors}

    if len(tumors) == 1:
        sys.stderr.write(log('# Intra-tumor analysis\n'))
        single(
            tumors[list(tumors)[0]],
            clones[list(clones)[0]],
            props[list(props)[0]],
            base[list(base)[0]],
            args,
        )
    else:
        sys.stderr.write(log('# Inter-tumors analysis\n'))
        multiple(tumors, clones, props, base, args)


def pp(tumor, clones, props, args):
    bases = {}
    for patient in tumor:
        sys.stderr.write(info('## PATIENT: {}\n'.format(patient)))
        counter = []
        for sample in props[patient]:
            purity = sum(
                float(props[patient][sample][i])
                for i in props[patient][sample]
                if i != 'normal'
            )
            scaled = {
                i: (props[patient][sample][i] / purity)
                if purity > 0.0
                else 0.0
                for i in props[patient][sample]
                if i != 'normal'
            }
            length = sum(
                float(s[1] - s[0])
                for c in tumor[patient]
                for s in tumor[patient][c]
            )
            ploidy = (
                sum(
                    float(sum(tumor[patient][c][s][i]))
                    * float(s[1] - s[0])
                    * scaled[i]
                    for c in tumor[patient]
                    for s in tumor[patient][c]
                    for i in scaled
                )
                / length
            )
            wgd = 2 if ploidy < args['threshold'] else 4
            sys.stderr.write(
                info(
                    '### SAMPLE: {} -- PURITY: {} -- PLOIDY: {} -- CLASSIFICATION: {}\n'.format(
                        sample,
                        purity,
                        ploidy,
                        'DIPLOID' if wgd == 2 else 'TETRAPLOID',
                    )
                )
            )
            counter.append(wgd)
        counter = Counter(counter)
        bases[patient] = argmax(counter)
    return bases


def single(tumor, clones, props, base, args):
    out = 'intratumor-clones-totalcn.pdf'
    sys.stderr.write(
        log(
            '# Plotting total copy-number clone profiles in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    profiles(tumor, clones, props, args, out)

    out = 'intratumor-clones-allelecn.pdf'
    sys.stderr.write(
        log(
            '# Plotting allele-specific copy-number clone profiles in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    allelicprofiles(tumor, clones, props, args, out)

    out = 'intratumor-copynumber-totalcn.pdf'
    sys.stderr.write(
        log(
            '# Plotting total copy-number proportions in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    cnproportions(tumor, base, clones, props, args, out)

    out = 'intratumor-copynumber-allelecn.pdf'
    sys.stderr.write(
        log(
            '# Plotting allele-specific copy-number proportions in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    allelicproportions(tumor, int(float(base) / 2), clones, props, args, out)

    out = 'intratumor-profiles.pdf'
    sys.stderr.write(
        log(
            '# Plotting clone profiles in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    gridprofiles(tumor, base, clones, props, args, out)

    out = 'intratumor-profilesreduced.pdf'
    sys.stderr.write(
        log(
            '# Plotting reduced-clone profiles in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    gridprofilesreduced(tumor, base, clones, props, args, out)

    # Run tumor mixture analyses if we have a multiple tumor samples (keys in props)
    if len(props) > 1:
        out = 'intratumor-mixtures.pdf'
        sys.stderr.write(
            log(
                '# Plotting reduced mixtures in {}\n'.format(
                    os.path.join(args['rundir'], out)
                )
            )
        )
        gridmixtures(tumor, base, clones, props, args, out)

        out = 'intratumor-subclonality.pdf'
        sys.stderr.write(
            log(
                '# Plotting reduced mixtures in {}\n'.format(
                    os.path.join(args['rundir'], out)
                )
            )
        )
        subclonal(tumor, base, clones, props, args, out)


def profiles(tumor, clones, props, args, out):
    proj = join(tumor, clones, args['profileres'])

    tclones = [i for i in clones if i != 'normal']
    shift = 0.1
    if len(tclones) % 2 == 0:
        available = [-shift / 2] + [
            -shift / 2 - shift * (i + 1)
            for i in range(int(len(tclones) / 2) - 1)
        ]
        available += [shift / 2] + [
            shift / 2 + shift * (i + 1)
            for i in range(int(len(tclones) / 2) - 1)
        ]
        available = iter(sorted(available))
        level = {clone: next(available) for clone in tclones}
    else:
        available = (
            [0]
            + [-shift * (i + 1) for i in range(int(len(tclones) / 2))]
            + [shift * (i + 1) for i in range(int(len(tclones) / 2))]
        )
        available = iter(sorted(available))
        level = {clone: next(available) for clone in tclones}
    pal = iter(
        map(
            mpl.colors.rgb2hex,
            args['clonepalette'](np.linspace(0, 1, len(tclones))),
        )
    )
    style = {name: next(pal) for i, name in enumerate(tclones)}

    plt.figure(figsize=args['clonefigsize'])
    pos = []
    x = 0
    for c in sorted(proj, key=sortchr):
        for s in sorted(proj[c], key=(lambda x: x[0])):
            pos.append((c, x))
            for i in proj[c][s]:
                if i in tclones:
                    y = sum(proj[c][s][i]) + level[i]
                    plt.scatter(x, y, c=style[i], marker='|', s=12)
            x += 1
    plt.xlim(xmin=0, xmax=x)
    ymin, ymax = plt.ylim()
    x = 0
    for c in sorted(proj, key=sortchr):
        for s in sorted(proj[c], key=(lambda x: x[0])):
            x += 1
        plt.plot((x, x), (0, ymax + 0.4), '--b', linewidth=0.2)
    addchrplt(pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')


def allelicprofiles(tumor, clones, props, args, out):
    proj = join(tumor, clones, args['profileres'])

    tclones = [i for i in clones if i != 'normal']
    shift = 0.1
    if len(tclones) % 2 == 0:
        available = [-shift / 2] + [
            -shift / 2 - shift * (i + 1)
            for i in range(int(len(tclones) / 2) - 1)
        ]
        available += [shift / 2] + [
            shift / 2 + shift * (i + 1)
            for i in range(int(len(tclones) / 2) - 1)
        ]
        available = iter(sorted(available))
        level = {clone: next(available) for clone in tclones}
    else:
        available = (
            [0]
            + [-shift * (i + 1) for i in range(int(len(tclones) / 2))]
            + [shift * (i + 1) for i in range(int(len(tclones) / 2))]
        )
        available = iter(sorted(available))
        level = {clone: next(available) for clone in tclones}
    pal = iter(
        map(
            mpl.colors.rgb2hex,
            args['clonepalette'](np.linspace(0, 1, len(tclones))),
        )
    )
    style = {name: next(pal) for i, name in enumerate(tclones)}

    plt.figure(figsize=args['clonefigsize'])

    x = 0
    pos = []
    for c in sorted(proj, key=sortchr):
        for s in sorted(proj[c], key=(lambda x: x[0])):
            pos.append((c, x))
            for i in proj[c][s]:
                if i != 'normal':
                    yA = proj[c][s][i][0] + level[i]
                    yB = -(proj[c][s][i][1] + level[i])
                    plt.scatter(x, yA, c=style[i], marker='|', s=12)
                    plt.scatter(x, yB, c=style[i], marker='|', s=12)
            x += 1

    plt.xlim(xmin=0, xmax=x)
    ymin, ymax = plt.ylim()
    ymin -= 0.4
    ymax += 0.4
    x = 0
    for c in sorted(proj, key=sortchr):
        for s in sorted(proj[c], key=(lambda x: x[0])):
            x += 1
        plt.plot((x, x), (ymin, ymax), '--b', linewidth=0.2)
    plt.fill_between(
        range(0, x),
        [0 for i in range(x)],
        [ymax for i in range(x)],
        facecolor='red',
        interpolate=True,
        alpha=0.05,
    )
    plt.fill_between(
        range(0, x),
        [ymin for i in range(x)],
        [0 for i in range(x)],
        facecolor='blue',
        interpolate=True,
        alpha=0.05,
    )
    addchrplt(pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')


def cnproportions(tumor, base, clones, props, args, out):
    with PdfPages(os.path.join(args['rundir'], out)) as pdf:
        for sample in props:
            sys.stderr.write(
                info('## Plotting for sample {}...\n'.format(sample))
            )
            proj = join(tumor, clones, args['cnres'])
            sumu = lambda cns, cn: sum(
                float(props[sample][i]) for i in cns if cns[i] == cn
            )
            merge = {
                c: {
                    s: {
                        sum(cn): sumu(proj[c][s], cn)
                        for cn in set(proj[c][s].values())
                    }
                    for s in proj[c]
                }
                for c in proj
            }

            pos = [
                (c, s)
                for c in sorted(merge, key=sortchr)
                for s in sorted(merge[c], key=(lambda x: x[0]))
            ]
            cns = sorted(
                set(cn for c in merge for s in merge[c] for cn in merge[c][s]),
                reverse=True,
            )
            pal2 = iter(sns.color_palette('coolwarm', 1))
            style = {base: next(pal2)}
            # palA = iter(sns.color_palette("Reds", len([cn for cn in cns if cn > base])))
            palA = iter(
                sns.color_palette(
                    'YlOrRd_r', len([cn for cn in cns if cn > base])
                )
            )
            for x in [cn for cn in cns if cn > base]:
                style[x] = next(palA)
            # palD = iter(sns.color_palette("Blues", len([cn for cn in cns if cn < base])))
            palD = iter(
                sns.color_palette(
                    'YlGnBu', len([cn for cn in cns if cn < base])
                )
            )
            for x in [cn for cn in cns if cn < base]:
                style[x] = next(palD)
            level = {s: 1.0 for s in pos}

            fig = plt.figure(figsize=args['propsfigsize'])

            for cn in cns:
                g = (
                    lambda s: merge[s[0]][s[1]][cn]
                    if cn in merge[s[0]][s[1]]
                    else 0.0
                )
                minv = lambda v: v if v > 0.02 else 0.0
                df = pd.DataFrame(
                    [
                        {
                            'Genome positions': x,
                            'Copy-number Proportions': minv(level[s]),
                        }
                        for x, s in enumerate(pos)
                    ]
                )
                sns.barplot(
                    data=df,
                    x='Genome positions',
                    y='Copy-number Proportions',
                    color=style[cn],
                    label=str(cn),
                )
                level = {s: level[s] - g(s) for s in pos}

            ticks = [
                (x, s[0])
                for x, s in enumerate(pos)
                if x == 0 or pos[x - 1][0] != pos[x][0]
            ]
            plt.xticks([x[0] for x in ticks], [x[1] for x in ticks])
            plt.legend(
                loc='center left',
                fancybox=True,
                shadow=True,
                bbox_to_anchor=(1, 0.5),
            )
            plt.ylim(ymin=0, ymax=1.0)
            fig.autofmt_xdate()
            pdf.savefig()
            plt.close()


def allelicproportions(tumor, base, clones, props, args, out):
    with PdfPages(os.path.join(args['rundir'], out)) as pdf:
        for sample in props:
            sys.stderr.write(
                info('## Plotting for sample {}...\n'.format(sample))
            )
            proj = join(
                tumor,
                [i for i in clones if props[sample][i] > 0.0],
                args['cnres'],
            )
            sumu = lambda cns, cn: sum(
                float(props[sample][i]) for i in cns if cns[i][0] == cn[0]
            )
            mergeA = {
                c: {
                    s: {
                        cn[0]: sumu(proj[c][s], cn)
                        for cn in set(proj[c][s].values())
                    }
                    for s in proj[c]
                }
                for c in proj
            }
            sumu = lambda cns, cn: sum(
                float(props[sample][i]) for i in cns if cns[i][1] == cn[1]
            )
            mergeB = {
                c: {
                    s: {
                        cn[1]: sumu(proj[c][s], cn)
                        for cn in set(proj[c][s].values())
                    }
                    for s in proj[c]
                }
                for c in proj
            }

            pos = [
                (c, s)
                for c in sorted(mergeA, key=sortchr)
                for s in sorted(mergeA[c], key=(lambda x: x[0]))
            ]
            cnsA = sorted(
                set(
                    cn
                    for c in mergeA
                    for s in mergeA[c]
                    for cn in mergeA[c][s]
                ),
                reverse=True,
            )
            cnsB = sorted(
                set(
                    cn
                    for c in mergeB
                    for s in mergeB[c]
                    for cn in mergeB[c][s]
                ),
                reverse=True,
            )
            cns = sorted(set(cnsA) | set(cnsB))

            pal2 = iter(sns.color_palette('coolwarm', 1))
            style = {base: next(pal2)}
            # palA = iter(sns.color_palette("Reds", len([cn for cn in cns if cn > base])))
            palA = iter(
                sns.color_palette(
                    'YlOrRd', len([cn for cn in cns if cn > base])
                )
            )
            for x in [cn for cn in cns if cn > base]:
                style[x] = next(palA)
            # palD = iter(sns.color_palette("Blues", len([cn for cn in cns if cn < base])))
            palD = iter(
                sns.color_palette(
                    'YlGnBu', len([cn for cn in cns if cn < base])
                )
            )
            for x in [cn for cn in cns if cn < base]:
                style[x] = next(palD)

            fig = plt.figure(figsize=args['propsfigsize'])

            # Plot allele A

            level = {s: 1.0 for s in pos}
            for cn in cns:
                g = (
                    lambda s: mergeA[s[0]][s[1]][cn]
                    if cn in mergeA[s[0]][s[1]]
                    else 0.0
                )
                minv = lambda v: v if v > 0.02 else 0.0
                df = pd.DataFrame(
                    [
                        {
                            'Genome positions': x,
                            'Copy-number Proportions': minv(level[s]),
                        }
                        for x, s in enumerate(pos)
                    ]
                )
                sns.barplot(
                    data=df,
                    x='Genome positions',
                    y='Copy-number Proportions',
                    color=style[cn],
                    label=str(cn),
                )
                level = {s: level[s] - g(s) for s in pos}

            # Plot allele B

            level = {s: -1.0 for s in pos}
            for cn in cns:
                g = (
                    lambda s: mergeB[s[0]][s[1]][cn]
                    if cn in mergeB[s[0]][s[1]]
                    else 0.0
                )
                minv = lambda v: v if v < 0.02 else 0.0
                df = pd.DataFrame(
                    [
                        {
                            'Genome positions': x,
                            'Copy-number Proportions': level[s],
                        }
                        for x, s in enumerate(pos)
                    ]
                )
                sns.barplot(
                    data=df,
                    x='Genome positions',
                    y='Copy-number Proportions',
                    color=style[cn],
                    label=str(cn),
                )
                level = {s: level[s] + g(s) for s in pos}

            ticks = [
                (x, s[0])
                for x, s in enumerate(pos)
                if x == 0 or pos[x - 1][0] != pos[x][0]
            ]
            plt.xticks([x[0] for x in ticks], [x[1] for x in ticks])
            plt.legend(
                loc='center left',
                fancybox=True,
                shadow=True,
                bbox_to_anchor=(1, 0.5),
            )
            plt.ylim(ymin=-1.0, ymax=1.0)
            fig.autofmt_xdate()
            pdf.savefig()
            plt.close()


def gridprofiles(tumor, base, clones, props, args, out):
    proj = join(tumor, clones, args['clusterres'])
    pos = [
        (c, s)
        for c in sorted(proj, key=sortchr)
        for s in sorted(proj[c], key=(lambda x: x[0]))
    ]
    pal_sample = sns.color_palette('YlGn', 10)
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {c: next(palette) for c in sorted(tumor, key=sortchr)}
    col_colors = {}
    row_colors = {}
    am = set()
    de = set()

    data = []
    for c in sorted([i for i in clones]):
        for x, s in enumerate(pos):
            cn = sum(proj[s[0]][s[1]][c])
            data.append({'Clone': c, 'Genome': x, 'Amp-Del': cn})
            if cn > base:
                am.add(cn)
            elif cn < base:
                de.add(cn)
            col_colors[x] = chr_colors[s[0]]
            row_colors[c] = {
                sam: pal_sample[min(9, int(round(props[sam][c] * 10)))]
                for sam in props
            }
    if len(am) == 0:
        am.add(base)
    if len(de) == 0:
        de.add(base)

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='Amp-Del',
        columns=['Genome'],
        index=['Clone'],
        aggfunc='first',
    )

    para = {}
    para['data'] = table
    para['cmap'] = 'coolwarm'
    para['center'] = base
    para['cbar_kws'] = {'ticks': range(min(de), max(am) + 1)}
    para['yticklabels'] = True
    para['row_cluster'] = True
    para['xticklabels'] = False
    para['col_cluster'] = False
    para['method'] = args['linkage']
    para['metric'] = cndistance
    para['figsize'] = args['clusterfigsize']
    para['col_colors'] = pd.DataFrame(
        [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
    ).set_index('index')
    para['row_colors'] = pd.DataFrame(
        [
            dict(list({'index': row}.items()) + list(row_colors[row].items()))
            for row in table.index
        ]
    ).set_index('index')
    g = sns.clustermap(**para)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def gridprofilesreduced(tumor, base, clones, props, args, out):
    proj = join(tumor, clones, args['clusterres'])
    red = reduction(proj, base)
    pos = [
        (c, s)
        for c in sorted(red, key=sortchr)
        for s in sorted(red[c], key=(lambda x: x[0]))
    ]
    pal_sample = sns.color_palette('YlGn', 10)
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {c: next(palette) for c in sorted(tumor, key=sortchr)}
    col_colors = {}
    row_colors = {}

    data = []
    for c in sorted([i for i in clones]):
        for x, s in enumerate(pos):
            data.append(
                {'Clone': c, 'Genome': x, 'Amp-Del': red[s[0]][s[1]][c]}
            )
            col_colors[x] = chr_colors[s[0]]
            row_colors[c] = {
                sam: pal_sample[min(9, int(round(props[sam][c] * 10)))]
                for sam in props
            }

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='Amp-Del',
        columns=['Genome'],
        index=['Clone'],
        aggfunc='first',
    )
    myColors = ('#67a9cf', '#f7f7f7', '#ef8a62')
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    para = {}
    para['data'] = table
    para['cmap'] = cmap
    para['cbar_kws'] = {
        'ticks': [-1, 0, 1],
        'boundaries': np.linspace(-1, 1, 4),
    }
    para['yticklabels'] = True
    para['row_cluster'] = True
    para['xticklabels'] = False
    para['col_cluster'] = False
    para['method'] = args['linkage']
    para['metric'] = cndistance
    para['figsize'] = args['clusterfigsize']
    para['col_colors'] = pd.DataFrame(
        [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
    ).set_index('index')
    para['row_colors'] = pd.DataFrame(
        [
            dict(list({'index': row}.items()) + list(row_colors[row].items()))
            for row in table.index
        ]
    ).set_index('index')
    g = sns.clustermap(**para)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def gridmixtures(tumor, base, clones, props, args, out):
    projp = join(tumor, clones, args['clusterres'])
    redp = reduction(projp, base)
    pos = [
        (c, s)
        for c in sorted(redp, key=sortchr)
        for s in sorted(redp[c], key=(lambda x: x[0]))
    ]
    pal_clone = sns.color_palette('YlGn', 10)
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {c: next(palette) for c in sorted(tumor, key=sortchr)}
    col_colors = {}
    row_colors = {}

    data = []
    for p in props:
        sumu = lambda cns, cn: sum(
            float(props[p][i]) for i in cns if cns[i] == cn
        )
        mergep = {
            c: {
                s: {
                    cn: sumu(redp[c][s], cn) for cn in set(redp[c][s].values())
                }
                for s in redp[c]
            }
            for c in redp
        }
        assert False not in set(
            0.99 <= sum(mergep[c][s][cn] for cn in mergep[c][s]) <= 1.01
            for c in mergep
            for s in mergep[c]
        )
        for x, s in enumerate(pos):
            value = sum(
                mergep[s[0]][s[1]][cn] * cn for cn in mergep[s[0]][s[1]]
            )
            data.append({'Sample': p, 'Genome': x, 'value': value})
            col_colors[x] = chr_colors[s[0]]
            row_colors[p] = {
                i: pal_clone[min(9, int(round(props[p][i] * 10)))]
                for i in clones
            }

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='value',
        columns=['Genome'],
        index=['Sample'],
        aggfunc='first',
    )

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        para = {}
        para['data'] = table
        para['cmap'] = 'coolwarm'  # "RdBu_r"
        # para['cbar_kws'] = {"ticks":[-2, -1, 0, 1, 2], "boundaries": np.linspace(-2, 2, 6)}
        para['yticklabels'] = True
        para['row_cluster'] = True
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['method'] = args['linkage']
        para['metric'] = similaritysample
        para['figsize'] = args['clusterfigsize']
        para['col_colors'] = pd.DataFrame(
            [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
        ).set_index('index')
        para['row_colors'] = pd.DataFrame(
            [
                dict(
                    list({'index': row}.items())
                    + list(row_colors[row].items())
                )
                for row in table.index
            ]
        ).set_index('index')
        g = sns.clustermap(**para)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def subclonal(tumor, base, clones, props, args, out):
    assert base in {2, 4}
    abase = 1 if base == 2 else 2
    proj = join(tumor, clones, args['clusterres'])
    pos = [
        (c, s)
        for c in sorted(proj, key=sortchr)
        for s in sorted(proj[c], key=(lambda x: x[0]))
    ]
    pal_clone = sns.color_palette('YlGn', 10)
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {c: next(palette) for c in sorted(tumor, key=sortchr)}
    col_colors = {}
    row_colors = {}

    data = []
    for p in props:
        for x, s in enumerate(pos):
            cns = set(
                proj[s[0]][s[1]][i] for i in proj[s[0]][s[1]] if i != 'normal'
            )
            merge = {
                cn: sum(
                    props[p][i]
                    for i in proj[s[0]][s[1]]
                    if proj[s[0]][s[1]][i] == cn and i != 'normal'
                )
                for cn in cns
            }
            cns = set(cn for cn in cns if merge[cn] >= args['minu'])
            if cns == {(abase, abase)}:
                value = 0
            elif False not in set(
                n[0] >= abase and n[1] >= abase for n in cns
            ):
                if len(cns) == 1:
                    value = 2
                else:
                    value = 1
            elif False not in set(
                n[0] <= abase and n[1] <= abase for n in cns
            ):
                if len(cns) == 1:
                    value = -2
                else:
                    value = -1
            else:
                if len(cns) == 1:
                    value = 4
                else:
                    value = 3
            data.append({'Sample': p, 'Genome': x, 'value': value})
            col_colors[x] = chr_colors[s[0]]
            row_colors[p] = {
                i: pal_clone[min(9, int(round(props[p][i] * 10)))]
                for i in clones
            }

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='value',
        columns=['Genome'],
        index=['Sample'],
        aggfunc='first',
    )
    myColors = (
        '#92c5de',
        '#0571b0',
        '#f7f7f7',
        '#ca0020',
        '#f4a582',
        '#7b3294',
        '#c2a5cf',
    )
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        para = {}
        para['data'] = table
        para['cmap'] = cmap
        labels = [
            'Clonal deletion',
            'Subclonal deletion',
            'Neutral',
            'Subclonal amplification',
            'Clonal amplification',
            'Sublonal mix',
            'Clonal mix',
        ]
        para['cbar_kws'] = {
            'ticks': [-2, -1, 0, 1, 2, 3, 4],
            'boundaries': np.linspace(-2, 4, 8),
        }
        para['yticklabels'] = True
        para['row_cluster'] = True
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['method'] = args['linkage']
        para['metric'] = similarity
        para['figsize'] = args['clusterfigsize']
        para['col_colors'] = pd.DataFrame(
            [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
        ).set_index('index')
        para['row_colors'] = pd.DataFrame(
            [
                dict(
                    list({'index': row}.items())
                    + list(row_colors[row].items())
                )
                for row in table.index
            ]
        ).set_index('index')
        g = sns.clustermap(**para)
        cax = plt.gcf().axes[-1]
        cax.set_yticklabels(labels)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def reduction(proj, base):
    classify = (
        lambda c, i: 0 if c == base or i == 'normal' else 1 if c > base else -1
    )
    reduction = {
        c: {
            s: {i: classify(sum(proj[c][s][i]), i) for i in proj[c][s]}
            for s in proj[c]
        }
        for c in proj
    }
    assert False not in set(
        reduction[c][s]['normal'] == 0 for c in proj for s in proj[c]
    )
    return reduction


def join(tumor, clones, resolution):
    proj = {}
    for c in tumor:
        bins = sorted(list(tumor[c]), key=(lambda x: x[0]))
        proj[c] = {}
        while bins:
            tmp = bins[:resolution]
            counts = {
                i: dict(Counter([tumor[c][s][i] for s in tmp])) for i in clones
            }
            proj[c][tmp[0][0], tmp[-1][1]] = {
                i: max(counts[i], key=(lambda x: counts[i][x])) for i in clones
            }
            bins = bins[resolution:]
    return proj


def multiple(tumor, clones, props, base, args):
    sys.stderr.write(log('# Uniforming bin segmention across all patients\n'))
    tumor = segmenting(tumor, clones, props)

    out = 'intertumors-profilesfull.pdf'
    sys.stderr.write(
        log(
            '# Plotting inter-tumors clone profiles in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    intergridfullprofiles(tumor, base, clones, props, args, out)

    out = 'intertumors-profilesreduced.pdf'
    sys.stderr.write(
        log(
            '# Plotting inter-tumors reduced-clone profiles in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    intergridreducedprofiles(tumor, base, clones, props, args, out)

    out = 'intertumors-mixtures.pdf'
    sys.stderr.write(
        log(
            '# Plotting inter-tumors mixtures in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    intergridsamplesclusters(tumor, base, clones, props, args, out)

    out = 'intertumors-subclonality.pdf'
    sys.stderr.write(
        log(
            '# Plotting inter-tumors subclonality in {}\n'.format(
                os.path.join(args['rundir'], out)
            )
        )
    )
    intergridsubclonality(tumor, base, clones, props, args, out)


def intergridfullprofiles(tumor, base, clones, props, args, out):
    proj = interjoin(tumor, clones, args['clusterres'])
    t = list(proj)[0]
    pos = [
        (c, s)
        for c in sorted(proj[t], key=sortchr)
        for s in sorted(proj[t][c], key=(lambda x: x[0]))
    ]
    palette = cycle(sns.color_palette('Pastel1', min(10, len(list(tumor)))))
    pat_colors = {pat: next(palette) for pat in tumor}
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {
        c: next(palette) for c in sorted(tumor[list(tumor)[0]], key=sortchr)
    }
    col_colors = {}
    row_colors = {}

    data = []
    for pat1 in proj:
        for c1 in [t for t in clones[pat1] if t != 'normal']:
            for x, s in enumerate(pos):
                data.append(
                    {
                        'Patient clone': '{}:{}'.format(pat1, c1),
                        'Genome': x,
                        'value': sum(proj[pat1][s[0]][s[1]][c1]),
                    }
                )
                col_colors[x] = chr_colors[s[0]]
                row_colors['{}:{}'.format(pat1, c1)] = pat_colors[pat1]

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='value',
        columns=['Genome'],
        index=['Patient clone'],
        aggfunc='first',
    )

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        para = {}
        para['data'] = table
        para['cmap'] = 'coolwarm'
        para['center'] = 2
        para['xticklabels'] = True
        para['yticklabels'] = True
        para['xticklabels'] = False
        para['row_cluster'] = False
        para['col_cluster'] = False
        para['cbar_kws'] = {'label': 'Total copy number'}
        para['figsize'] = args['clusterfigsize']
        para['method'] = args['linkage']
        para['metric'] = cndistance
        para['col_colors'] = pd.DataFrame(
            [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
        ).set_index('index')
        para['row_colors'] = pd.DataFrame(
            [{'index': row, 'patient': row_colors[row]} for row in table.index]
        ).set_index('index')
        g = sns.clustermap(**para)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def intergridreducedprofiles(tumor, base, clones, props, args, out):
    proj = interjoin(tumor, clones, args['clusterres'])
    red = interreduction(proj, base)
    t = list(red)[0]
    pos = [
        (c, s)
        for c in sorted(red[t], key=sortchr)
        for s in sorted(red[t][c], key=(lambda x: x[0]))
    ]
    palette = cycle(sns.color_palette('Pastel1', min(9, len(list(tumor)))))
    pat_colors = {pat: next(palette) for pat in tumor}
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {
        c: next(palette) for c in sorted(tumor[list(tumor)[0]], key=sortchr)
    }
    col_colors = {}
    row_colors = {}

    data = []
    for pat1 in red:
        for c1 in [t for t in clones[pat1] if t != 'normal']:
            for x, s in enumerate(pos):
                data.append(
                    {
                        'Patient clone': '{}:{}'.format(pat1, c1),
                        'Genome': x,
                        'value': red[pat1][s[0]][s[1]][c1],
                    }
                )
                col_colors[x] = chr_colors[s[0]]
                row_colors['{}:{}'.format(pat1, c1)] = pat_colors[pat1]

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='value',
        columns=['Genome'],
        index=['Patient clone'],
        aggfunc='first',
    )

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        para = {}
        para['data'] = table
        para['cmap'] = 'coolwarm'
        para['xticklabels'] = True
        para['yticklabels'] = True
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['figsize'] = args['clusterfigsize']
        para['cbar_kws'] = {
            'ticks': [-1, 0, 1],
            'boundaries': np.linspace(-1, 1, 4),
        }
        para['method'] = args['linkage']
        para['metric'] = similarity
        para['col_colors'] = pd.DataFrame(
            [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
        ).set_index('index')
        para['row_colors'] = pd.DataFrame(
            [{'index': row, 'patient': row_colors[row]} for row in table.index]
        ).set_index('index')
        g = sns.clustermap(**para)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def intergridsamplesclusters(tumor, base, clones, props, args, out):
    data = []
    proj = interjoin(tumor, clones, args['clusterres'])
    red = interreduction(proj, base)
    t = list(red)[0]
    pos = [
        (c, s)
        for c in sorted(red[t], key=sortchr)
        for s in sorted(red[t][c], key=(lambda x: x[0]))
    ]
    palette = cycle(sns.color_palette('Pastel1', min(9, len(list(tumor)))))
    pat_colors = {pat: next(palette) for pat in tumor}
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {
        c: next(palette) for c in sorted(tumor[list(tumor)[0]], key=sortchr)
    }
    col_colors = {}
    row_colors = {}

    for pat in tumor:
        for p in props[pat]:
            sumu = lambda cns, cn: sum(
                float(props[pat][p][i]) if i in props[pat][p] else 0.0
                for i in cns
                if cns[i] == cn
            )
            mergep = {
                c: {
                    s: {
                        cn: sumu(red[pat][c][s], cn)
                        for cn in set(red[pat][c][s].values())
                    }
                    for s in red[pat][c]
                }
                for c in red[pat]
            }
            assert False not in set(
                0.99 <= sum(mergep[c][s][cn] for cn in mergep[c][s]) <= 1.01
                for c in mergep
                for s in mergep[c]
            )

            for x, s in enumerate(pos):
                value = sum(
                    mergep[s[0]][s[1]][cn] * cn for cn in mergep[s[0]][s[1]]
                )
                data.append(
                    {
                        'Patient sample': '{}-{}'.format(pat, p),
                        'Genome': x,
                        'value': value,
                    }
                )
                col_colors[x] = chr_colors[s[0]]
                row_colors['{}-{}'.format(pat, p)] = pat_colors[pat]

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='value',
        columns=['Genome'],
        index=['Patient sample'],
        aggfunc='first',
    )

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        para = {}
        para['data'] = table
        para['cmap'] = 'coolwarm'
        para['xticklabels'] = True
        para['yticklabels'] = True
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['figsize'] = args['clusterfigsize']
        para['method'] = args['linkage']
        para['metric'] = similaritysample
        para['col_colors'] = pd.DataFrame(
            [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
        ).set_index('index')
        para['row_colors'] = pd.DataFrame(
            [{'index': row, 'patient': row_colors[row]} for row in table.index]
        ).set_index('index')
        g = sns.clustermap(**para)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def intergridsubclonality(tumor, base, clones, props, args, out):
    data = []
    proj = interjoin(tumor, clones, args['clusterres'])
    t = list(proj)[0]
    pos = [
        (c, s)
        for c in sorted(proj[t], key=sortchr)
        for s in sorted(proj[t][c], key=(lambda x: x[0]))
    ]
    palette = cycle(sns.color_palette('Pastel1', min(9, len(list(tumor)))))
    pat_colors = {pat: next(palette) for pat in tumor}
    palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {
        c: next(palette) for c in sorted(tumor[list(tumor)[0]], key=sortchr)
    }
    col_colors = {}
    row_colors = {}

    for pat in tumor:
        assert base[pat] == 2 or base[pat] == 4
        abase = 1 if base[pat] == 2 else 2
        for p in props[pat]:
            for x, s in enumerate(pos):
                cns = set(
                    proj[pat][s[0]][s[1]][i]
                    for i in proj[pat][s[0]][s[1]]
                    if i != 'normal'
                )
                merge = {
                    cn: sum(
                        props[pat][p][i]
                        for i in proj[pat][s[0]][s[1]]
                        if proj[pat][s[0]][s[1]][i] == cn and i != 'normal'
                    )
                    for cn in cns
                }
                cns = set(cn for cn in cns if merge[cn] >= args['minu'])
                if cns == {(abase, abase)}:
                    value = 0
                elif False not in set(
                    n[0] >= abase and n[1] >= abase for n in cns
                ):
                    if len(cns) == 1:
                        value = 2
                    else:
                        value = 1
                elif False not in set(
                    n[0] <= abase and n[1] <= abase for n in cns
                ):
                    if len(cns) == 1:
                        value = -2
                    else:
                        value = -1
                else:
                    if len(cns) == 1:
                        value = 4
                    else:
                        value = 3
                data.append(
                    {
                        'Patient sample': '{}-{}'.format(pat, p),
                        'Genome': x,
                        'value': value,
                    }
                )
                col_colors[x] = chr_colors[s[0]]
                row_colors['{}-{}'.format(pat, p)] = pat_colors[pat]

    df = pd.DataFrame(data)
    table = pd.pivot_table(
        df,
        values='value',
        columns=['Genome'],
        index=['Patient sample'],
        aggfunc='first',
    )
    myColors = (
        '#92c5de',
        '#0571b0',
        '#f7f7f7',
        '#ca0020',
        '#f4a582',
        '#7b3294',
        '#c2a5cf',
    )
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        para = {}
        para['data'] = table
        para['cmap'] = cmap
        labels = [
            'Clonal deletion',
            'Subclonal deletion',
            'Neutral',
            'Subclonal amplification',
            'Clonal amplification',
            'Sublonal mix',
            'Clonal mix',
        ]
        para['cbar_kws'] = {
            'ticks': [-2, -1, 0, 1, 2, 3, 4],
            'boundaries': np.linspace(-2, 4, 8),
        }
        para['xticklabels'] = True
        para['yticklabels'] = True
        para['xticklabels'] = False
        para['col_cluster'] = False
        para['figsize'] = args['clusterfigsize']
        para['method'] = args['linkage']
        para['metric'] = similaritysample
        para['col_colors'] = pd.DataFrame(
            [{'index': s, 'chromosomes': col_colors[s]} for s in table.columns]
        ).set_index('index')
        para['row_colors'] = pd.DataFrame(
            [{'index': row, 'patient': row_colors[row]} for row in table.index]
        ).set_index('index')
        g = sns.clustermap(**para)
        cax = plt.gcf().axes[-1]
        cax.set_yticklabels(labels)

    addchr(g, pos)
    plt.savefig(os.path.join(args['rundir'], out), bbox_inches='tight')
    plt.close()


def segmenting(tumor, clones, props):
    numpat = len(tumor)
    cbk = lambda c: set(b for pat in tumor for s in tumor[pat][c] for b in s)
    bk = {c: cbk(c) for c in set(c for pat in tumor for c in tumor[pat])}
    bk = {c: sorted(bk[c]) for c in bk}
    counts = {c: {b: 0 for b in zip(bk[c][:-1], bk[c][1:])} for c in bk}
    maps = {
        pat: {c: {b: None for b in counts[c]} for c in counts} for pat in tumor
    }

    def select(d, m, counts, breakpoints, numpat):
        for c in d:
            bk = deque(breakpoints[c])
            left = -1
            right = bk.popleft()
            for (_left, _right) in sorted(d[c], key=(lambda x: x[0])):
                while right != _right:
                    left = right
                    right = bk.popleft()
                    if _left <= left and right <= _right:
                        counts[c][left, right] += 1
                        assert counts[c][left, right] <= numpat
                        m[c][left, right] = (_left, _right)

    for pat in tumor:
        select(tumor[pat], maps[pat], counts, bk, numpat)

    taken = {
        c: set(b for b in counts[c] if counts[c][b] == numpat) for c in counts
    }
    tottaken = sum(
        1.0 for c in counts for b in counts[c] if counts[c][b] == numpat
    )
    tot = sum(1.0 for c in counts for b in counts[c])
    sys.stderr.write(
        info(
            '## Proportion of common bins kept: {}%\n'.format(
                (tottaken / tot * 100) if tot > 0 else 1
            )
        )
    )

    return {
        pat: {
            c: {b: tumor[pat][c][maps[pat][c][b]] for b in taken[c]}
            for c in taken
        }
        for pat in tumor
    }


def interreduction(proj, base):
    classify = (
        lambda c, pat, i: 0
        if c == base[pat] or i == 'normal'
        else 1
        if c > base[pat]
        else -1
    )
    reduction = {
        pat: {
            c: {
                s: {
                    i: classify(sum(proj[pat][c][s][i]), pat, i)
                    for i in proj[pat][c][s]
                }
                for s in proj[pat][c]
            }
            for c in proj[pat]
        }
        for pat in proj
    }
    assert False not in set(
        reduction[pat][c][s]['normal'] == 0
        for pat in proj
        for c in proj[pat]
        for s in proj[pat][c]
    )
    return reduction


def interjoin(tumor, clones, resolution):
    proj = {}
    for pat in tumor:
        proj[pat] = {}
        for c in tumor[pat]:
            bins = sorted(list(tumor[pat][c]), key=(lambda x: x[0]))
            proj[pat][c] = {}
            while bins:
                tmp = bins[:resolution]
                counts = {
                    i: dict(Counter([tumor[pat][c][s][i] for s in tmp]))
                    for i in clones[pat]
                }
                proj[pat][c][tmp[0][0], tmp[-1][1]] = {
                    i: max(counts[i], key=(lambda x: counts[i][x]))
                    for i in clones[pat]
                }
                bins = bins[resolution:]
    return proj


def cndistance(u, v):
    diff = list(u - v)
    amps = [abs(x) if x > 0 else 0 for x in diff]
    dels = [abs(x) if x < 0 else 0 for x in diff]
    dist = sum(max(amps[i] - amps[i - 1], 0) for i, x in enumerate(amps))
    dist += sum(max(dels[i] - dels[i - 1], 0) for i, x in enumerate(dels))
    return dist


def similarity(u, v):
    a = float(
        sum(u[i] == v[i] and u[i] != 0.0 and v[i] != 0 for i in range(len(u)))
    )
    b = float(sum(u[i] != 0 or v[i] != 0 for i in range(len(u))))
    return (a / b) if b > 0 else 0


def similaritysample(u, v):
    bothamp = lambda x, y: x > 0.0 and y > 0.0
    bothdel = lambda x, y: x < 0.0 and y < 0.0
    a = float(
        sum(
            (bothamp(u[i], v[i]) or bothdel(u[i], v[i]))
            and u[i] != 0
            and v[i] != 0
            for i in range(len(u))
        )
    )
    b = float(sum(u[i] != 0 or v[i] != 0 for i in range(len(u))))
    return (a / b) if b > 0 else 0


def readUCN(inputs, patnames):
    tumors = {}
    samples = {}
    clones = {}
    props = {}
    for fil in inputs:
        sys.stderr.write(
            info('## Reading {} as {}...\n'.format(fil, patnames[fil]))
        )
        with open(fil, 'r') as f:
            patient = patnames[fil]
            tumors[patient] = {}
            samples[patient] = set()
            clones[patient] = []
            props[patient] = {}
            for line in f:
                if len(clones[patient]) == 0:
                    assert line[0] == '#'
                    clones[patient] = [
                        f.split('_')[1]
                        for i, f in enumerate(line.strip().split()[11:])
                        if i % 2 == 0
                    ]
                    if 'normal' not in clones[patient]:
                        raise ValueError(
                            error(
                                'normal is not present as a clone in {}'.format(
                                    patient
                                )
                            )
                        )
                else:
                    if len(line) > 1 and line[0] != '#':
                        parsed = line.strip().split()
                        sample = parsed[3]
                        samples[patient].add(sample)
                        chro = parsed[0]
                        if chro not in tumors[patient]:
                            tumors[patient][chro] = {}
                        start = int(parsed[1])
                        end = int(parsed[2])
                        pair = lambda l: [
                            (tuple(map(int, e.split('|'))), float(l[i + 1]))
                            for i, e in enumerate(l)
                            if i % 2 == 0
                        ]
                        read = {
                            clones[patient][i]: e[0]
                            for i, e in enumerate(pair(parsed[11:]))
                        }
                        if (start, end) not in tumors[patient][chro]:
                            tumors[patient][chro][start, end] = read
                        else:
                            for i in read:
                                read[i] == tumors[patient][chro][start, end][i]
                        check = {
                            clones[patient][i]: e[1]
                            for i, e in enumerate(pair(parsed[11:]))
                        }
                        if sample in props[patient]:
                            for i in clones[patient]:
                                assert check[i] == props[patient][sample][i]
                        else:
                            props[patient][sample] = {
                                i: check[i] for i in clones[patient]
                            }
                        assert (
                            0.999
                            <= sum(check[i] for i in clones[patient])
                            <= 1.001
                        )

    return tumors, clones, props


def addchrplt(pos):
    corners = []
    prev = 0
    val = pos[0][0]
    for x, s in enumerate(pos):
        if x != 0 and pos[x - 1][0] != pos[x][0]:
            corners.append((prev, x, val))
            prev = x
            val = s[0]
    corners.append((prev, x, val))
    ticks = [(int(float(o[1] + o[0] + 1) / 2.0), o[2]) for o in corners]
    plt.xticks(
        [x[0] for x in ticks], [x[1] for x in ticks], rotation=45, ha='center'
    )
    plt.yticks(rotation=0)


def addchr(g, pos, color=None):
    corners = []
    prev = 0
    for x, s in enumerate(pos):
        if x != 0 and pos[x - 1][0] != pos[x][0]:
            corners.append((prev, x))
            prev = x
    corners.append((prev, x))
    ax = g.ax_heatmap
    ticks = []
    for o in corners:
        ax.set_xticks(
            np.append(ax.get_xticks(), int(float(o[1] + o[0] + 1) / 2.0))
        )
        ticks.append(pos[o[0]][0])
    ax.set_xticklabels(ticks, rotation=45, ha='center')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)


def sortchr(x):
    if x.endswith('X'):
        return 23
    elif x.endswith('Y'):
        return 24
    else:
        return int(''.join([d for d in x if d.isdigit()]))


def argmax(d):
    return max(d, key=(lambda x: d[x]))


def argmin(d):
    return min(d, key=(lambda x: d[x]))


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def error(msg):
    return '{}{}{}'.format('\033[91m\033[1m', msg, '\033[0m')


def warning(msg):
    return '{}{}{}'.format('\033[93m\033[1m', msg, '\033[0m')


def log(msg):
    return '{}{}{}'.format('\033[95m\033[1m', msg, '\033[0m')


def info(msg):
    return '{}{}{}'.format('\033[96m', msg, '\033[0m')


def debug(msg):
    return '{}{}{}'.format('\033[92m', msg, '\033[0m')


if __name__ == '__main__':
    main()
