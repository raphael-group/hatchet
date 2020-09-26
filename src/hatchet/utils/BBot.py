#!/usr/bin/python2

import matplotlib as mpl
mpl.use('Agg')
import sys, os, argparse
import math
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
from scipy.stats import beta
import matplotlib.ticker as ticker
from scipy.stats import gaussian_kde
import matplotlib.colors as col
from matplotlib.pyplot import cm
from itertools import cycle
from collections import Counter

from ArgParsing import parse_bbot_args

plt.style.use('ggplot')
sns.set_style("whitegrid")


def main(args=None):
    sys.stderr.write(log("# Parsing and checking input arguments\n"))
    args = parse_bbot_args(args)
    sys.stdout.write(info("\n".join(["## {}:\t{}".format(key, args[key]) for key in args]) + '\n'))

    sys.stderr.write(log("# Reading input BBC file\n"))
    bbc, clusters = readBBC(args['input'])

    if args['fontscale'] != 1:
        sns.set(font_scale = args['fontscale'])

    if args['resolution'] is not None:
        sys.stderr.write(log("# Merging bins according to resolution\n"))
        bbc, clusters = join(bbc, clusters, args['resolution'])

    if args['ct'] is not None or args['st'] is not None:
        sys.stderr.write(log("# Bin's clusters are selected accordingly to the provided thresholds\n"))
        bbc, clusters = select(bbc, clusters, args)

    if args['command'] is None or args['command'] == 'RD':
        out = os.path.join(args['x'], 'readdepthratio.pdf')
        sys.stderr.write(log("# [RD] Plotting read-depth ratio (RDR) for all samples in {}\n".format(out)))
        rdr(bbc, args, out)

    if args['command'] is None or args['command'] == 'CRD':
        out = os.path.join(args['x'], 'readdepthratio_clustered.pdf')
        sys.stderr.write(log("# [CRD] Plotting the clustered read-depth ratio (RDR) for each sample in {}\n".format(out)))
        clurdr(bbc, clusters, args, out)

    if args['command'] is None or args['command'] == 'BAF':
        out = os.path.join(args['x'], 'ballelefrequency.pdf')
        sys.stderr.write(log("# [BAF] Plotting B-allele frequency (BAF) for all samples in {}\n".format(out)))
        baf(bbc, args, out)

    if args['command'] is None or args['command'] == 'CBAF':
        out = os.path.join(args['x'], 'ballelefrequency_clustered.pdf')
        sys.stderr.write(log("# [CBAF] Plotting the clustered B-allele frequency (BAF) for each sample in {}\n".format(out)))
        clubaf(bbc, clusters, args, out)

    if args['command'] is None or args['command'] == 'BB':
        out = os.path.join(args['x'], 'bb.pdf')
        sys.stderr.write(log("# [BB] Plotting RDR-BB for all samples in {}\n".format(out)))
        bb(bbc, clusters, args, out)

    if args['command'] is None or args['command'] == 'CBB':
        out = os.path.join(args['x'], 'bb_clustered.pdf')
        sys.stderr.write(log("# [CBB] Plotting clustered RDR-BB for all samples in {}\n".format(out)))
        clubb(bbc, clusters, args, out)

    if args['command'] is None or args['command'] == 'CLUSTER':
        if args['segfile'] is not None:
            seg = readSEG(args['segfile'])
            out = os.path.join(args['x'], 'clusters.pdf')
            sys.stderr.write(log("# [CLUSTER] Plotting clusters for all samples in {}\n".format(out)))
            clus(seg, args, out)
        else:
            sys.stderr.write(warning('### Provide a .seg file to also plot CLUSTER\n'))


def rdr(bbc, args, out):
    pos = [(c, s) for c in sorted(bbc, key=sortchr) for s in sorted(bbc[c], key=(lambda z : z[0]))]
    lx = 'Genome'
    ly = 'Read-depth ratio (RDR)'
    lh = 'Sample'
    data = [{lx : x, ly : bbc[b[0]][b[1]][p]['RDR'], lh : p} for x, b in enumerate(pos) for p in bbc[b[0]][b[1]]]
    df = pd.DataFrame(data)
    df.sort_values([lx, lh], ascending=[True, False])
    figsize = args['figsize'] if args['figsize'] is not None else (8, 2)
    s = args['markersize'] if args['markersize'] > 0 else 20
    mpl.rcParams['figure.figsize'] = (figsize[0], figsize[1])

    with PdfPages(out) as pdf:
        sys.stderr.write(info("## Plotting for all samples..\n"))
        g = sns.lmplot(data=df, x=lx, y=ly, hue=lh, palette=args['cmap'], fit_reg=False, size=figsize[0], aspect=figsize[1], markers="|", scatter_kws={"s":s})
        addchr(pos)
        coordinates(args, g)
        plt.xlim(xmin=0, xmax=(len(pos)+1))
        coordinates(args, g)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        for sample, group in df.groupby(lh):
            sys.stderr.write(info("## Plotting for {}..\n".format(sample)))
            g = sns.lmplot(data=group, x=lx, y=ly, fit_reg=False, size=figsize[0], aspect=figsize[1], markers="|", scatter_kws={"s":s})
            addchr(pos)
            coordinates(args, g)
            plt.title("Read-depth ratio in {}".format(sample))
            plt.xlim(xmin=0, xmax=(len(pos)+1))
            coordinates(args, g)
            pdf.savefig(bbox_inches='tight')
            plt.close()


def clurdr(bbc, clusters, args, out):
    pos = [(c, s) for c in sorted(bbc, key=sortchr) for s in sorted(bbc[c], key=(lambda z : z[0]))]
    lx = 'Genome'
    ly = 'Read-depth ratio (RDR)'
    g = 'Sample'
    lh = 'Cluster'
    data = [{lx : x, ly : bbc[b[0]][b[1]][p]['RDR'], g : p, lh : clusters[b[0]][b[1]]} for x, b in enumerate(pos) for p in bbc[b[0]][b[1]]]
    df = pd.DataFrame(data)
    df.sort_values([lx, lh], ascending=[True, True])
    figsize = args['figsize'] if args['figsize'] is not None else (8, 2)
    s = args['markersize'] if args['markersize'] > 0 else 20
    mpl.rcParams['figure.figsize'] = (figsize[0], figsize[1])

    with PdfPages(out) as pdf:
        for sample, group in df.groupby(g):
            sys.stderr.write(info("## Plotting for {}..\n".format(sample)))
            sns.lmplot(data=group, x=lx, y=ly, hue=lh, fit_reg=False, size=figsize[0], aspect=figsize[1], markers="|", scatter_kws={"s":s}, legend=False)
            addchr(pos)
            coordinates(args)
            plt.title("Read-depth ratio in {}".format(sample))
            plt.xlim(xmin=0, xmax=(len(pos)+1))
            pdf.savefig(bbox_inches='tight')
            plt.close()


def baf(bbc, args, out):
    pos = [(c, s) for c in sorted(bbc, key=sortchr) for s in sorted(bbc[c], key=(lambda z : z[0]))]
    lx = 'Genome'
    ly = 'B-allele frequency (BAF)'
    lh = 'Sample'
    data = [{lx : x, ly : bbc[b[0]][b[1]][p]['BAF'], lh : p} for x, b in enumerate(pos) for p in bbc[b[0]][b[1]]]
    df = pd.DataFrame(data)
    df.sort_values([lx, lh], ascending=[True, False])
    figsize = args['figsize'] if args['figsize'] is not None else (8, 2)
    s = args['markersize'] if args['markersize'] > 0 else 20
    mpl.rcParams['figure.figsize'] = (figsize[0], figsize[1])

    with PdfPages(out) as pdf:
        sys.stderr.write(info("## Plotting for all samples..\n"))
        g = sns.lmplot(data=df, x=lx, y=ly, hue=lh, palette=args['cmap'], fit_reg=False, size=figsize[0], aspect=figsize[1], markers="|", scatter_kws={"s":s})
        plt.ylim(ymax=0.5)
        addchr(pos)
        coordinates(args)
        plt.xlim(xmin=0, xmax=(len(pos)+1))
        pdf.savefig(bbox_inches='tight')
        plt.close()

        for sample, group in df.groupby(lh):
            sys.stderr.write(info("## Plotting for {}..\n".format(sample)))
            sns.lmplot(data=group, x=lx, y=ly, fit_reg=False, size=figsize[0], aspect=figsize[1], markers="|", scatter_kws={"s":s})
            plt.ylim(ymax=0.5)
            addchr(pos)
            coordinates(args)
            plt.title("B-allele frequency in {}".format(sample))
            plt.xlim(xmin=0, xmax=(len(pos)+1))
            pdf.savefig(bbox_inches='tight')
            plt.close()


def clubaf(bbc, clusters, args, out):
    pos = [(c, s) for c in sorted(bbc, key=sortchr) for s in sorted(bbc[c], key=(lambda z : z[0]))]
    lx = 'Genome'
    ly = 'B-allele frequency (BAF)'
    g = 'Sample'
    lh = 'Cluster'
    data = [{lx : x, ly : bbc[b[0]][b[1]][p]['BAF'], g : p, lh : clusters[b[0]][b[1]]} for x, b in enumerate(pos) for p in bbc[b[0]][b[1]]]
    df = pd.DataFrame(data)
    df.sort_values([lx, lh], ascending=[True, True])
    figsize = args['figsize'] if args['figsize'] is not None else (8, 2)
    s = args['markersize'] if args['markersize'] > 0 else 20
    mpl.rcParams['figure.figsize'] = (figsize[0], figsize[1])

    with PdfPages(out) as pdf:
        for sample, group in df.groupby(g):
            sys.stderr.write(info("## Plotting for {}..\n".format(sample)))
            sns.lmplot(data=group, x=lx, y=ly, hue=lh, fit_reg=False, size=figsize[0], aspect=figsize[1], markers="|", scatter_kws={"s":s}, legend=False)
            plt.ylim(ymax=0.5)
            addchr(pos)
            coordinates(args)
            plt.title("B-allele frequency in {}".format(sample))
            plt.xlim(xmin=0, xmax=(len(pos)+1))
            pdf.savefig(bbox_inches='tight')
            plt.close()


def bb(bbc, clusters, args, out):
    pos = [(c, s) for c in sorted(bbc, key=sortchr) for s in sorted(bbc[c], key=(lambda z : z[0]))]
    lx = 'RDR'
    ly = '0.5 - BAF'
    g = 'Sample'
    data = [{lx : bbc[b[0]][b[1]][p]['RDR'], ly : 0.5 - bbc[b[0]][b[1]][p]['BAF'], g : p} for b in pos for p in bbc[b[0]][b[1]]]
    df = pd.DataFrame(data)
    figsize = args['figsize'] if args['figsize'] is not None else (16, 8)
    s = args['markersize'] if args['markersize'] > 0 else 10

    with PdfPages(out) as pdf:
        for sample, group in df.groupby(g):
            sys.stderr.write(info("## Plotting for {}..\n".format(sample)))
            rdratio = np.array(group[lx])
            baf = np.array(group[ly])

            # Calculate the point density
            xy = np.vstack([rdratio,baf])
            z = gaussian_kde(xy)(xy)

            # Sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            rdratio, baf, z = rdratio[idx], baf[idx], z[idx]

            fig, ax = plt.subplots(1, figsize=figsize)
            cax = ax.scatter(rdratio, baf, c=z, cmap=plt.cm.jet, norm=col.LogNorm(),s=s)
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
            ax.grid(True)
            plt.colorbar(cax)
            plt.title(sample)
            coordinates(args)
            pdf.savefig(bbox_inches='tight')
            plt.close()


def clubb(bbc, clusters, args, out):
    pos = [(c, s) for c in sorted(bbc, key=sortchr) for s in sorted(bbc[c], key=(lambda z : z[0]))]
    lx = 'RDR'
    ly = '0.5 - BAF'
    g = 'Sample'
    lh = 'Cluster'
    size = {i : float(sum(clusters[b[0]][b[1]] == i for b in pos)) for i in set(clusters[b[0]][b[1]] for b in pos)}
    data = [{lx : bbc[b[0]][b[1]][p]['RDR'], ly : 0.5 - bbc[b[0]][b[1]][p]['BAF'], g : p, lh : clusters[b[0]][b[1]], 'size' : size[clusters[b[0]][b[1]]]} for b in pos for p in bbc[b[0]][b[1]]]
    df = pd.DataFrame(data)
    order = sorted(set(df[lh]), key=(lambda x : size[x]), reverse=True)
    figsize = args['figsize'] if args['figsize'] is not None else (8, 1.5)
    s = args['markersize'] if args['markersize'] > 0 else 20

    #with PdfPages(out) as pdf:
    #    for sample, group in df.groupby(g):
    #sys.stderr.write(info("## Plotting for {}..\n".format(sample)))
    if args['colwrap'] > 1:
        g = sns.lmplot(data=df, x=lx, y=ly, hue=lh, hue_order=order, palette=args['cmap'], fit_reg=False, size=figsize[0], aspect=figsize[1], scatter_kws={"s":s}, legend=False, col=g, col_wrap=args['colwrap'])
    else:
        g = sns.lmplot(data=df, x=lx, y=ly, hue=lh, hue_order=order, palette=args['cmap'], fit_reg=False, size=figsize[0], aspect=figsize[1], scatter_kws={"s":s}, legend=False, row=g)
    #plt.title("{}".format(sample))
    coordinates(args, g)
    #pdf.savefig(bbox_inches='tight')
    plt.savefig(out, bbox_inches='tight')
    plt.close()


def clus(seg, args, out):
    lx = 'Read-depth ratio (RDR)'
    ly = '0.5 - B-allele frequency (BAF)'
    g = 'Sample'
    lh = 'Cluster'
    samples = set(seg[list(seg)[0]])
    figsize = args['figsize'] if args['figsize'] is not None else (16, 10)
    s = args['markersize'] if args['markersize'] > 0 else 20
    mpl.rcParams['figure.figsize'] = (figsize[0], figsize[1])
    pal = cycle(sns.color_palette(args['cmap'], min(20, len(set(seg)))))
    col = {idx : next(pal) for idx in seg}

    with PdfPages(out) as pdf:
        for p in samples:
            sys.stderr.write(info("## Plotting for {}..\n".format(p)))
            for idx in seg:
                plt.scatter(seg[idx][p]['RDR'], 0.5 - seg[idx][p]['BAF'], c=col[idx], s=(seg[idx][p]['SIZE']**0.5)*20, alpha=0.8)
            plt.title("{}".format(p))
            coordinates(args)
            pdf.savefig(bbox_inches='tight')
            plt.close()


def readBBC(inp):
    bbc = {}
    clusters = {}
    samples = set()
    with open(inp, 'r') as f:
        for line in f:
            if len(line) > 0 and line[0] != '#':
                parsed = line.strip().split()
                chro = parsed[0]
                if chro not in bbc:
                    bbc[chro] = {}
                    clusters[chro] = {}
                start = int(parsed[1])
                end = int(parsed[2])
                if (start, end) not in bbc[chro]:
                    bbc[chro][start, end] = {}
                    clusters[chro][start, end] = {}
                clusters[chro][start, end] = parsed[10]
                sample = parsed[3]
                samples.add(sample)
                if sample not in bbc[chro][start, end]:
                    bbc[chro][start, end][sample] = {}
                bbc[chro][start, end][sample]['RDR'] = float(parsed[4])
                bbc[chro][start, end][sample]['BAF'] = float(parsed[9])
                bbc[chro][start, end][sample]['SNPS'] = float(parsed[5])

    for chro in bbc:
        for s in bbc[chro]:
            assert len(bbc[chro][s]) == len(samples)

    return bbc, clusters


def readSEG(inp):
    seg = {}
    samples = set()
    with open(inp, 'r') as f:
        for line in f:
            if len(line) > 0 and line[0] != '#':
                parsed = line.strip().split()
                idx = parsed[0]
                if idx not in seg:
                    seg[idx] = {}
                sample = parsed[1]
                samples.add(sample)
                if sample not in seg[idx]:
                    seg[idx][sample] = {}
                seg[idx][sample]['SIZE'] = int(parsed[2])
                seg[idx][sample]['RDR'] = float(parsed[3])
                seg[idx][sample]['BAF'] = float(parsed[8])

    for idx in seg:
        assert len(seg[idx]) == len(samples)

    return seg


def join(bbc, clusters, resolution):
    projbbc = {}
    projclu = {}
    samples = set(p for c in bbc for s in bbc[c] for p in bbc[c][s])
    for c in bbc:
        bins = sorted(list(bbc[c]), key=(lambda x : x[0]))
        projbbc[c] = {}
        projclu[c] = {}
        while bins:
            tmp = bins[:resolution]
            projbbc[c][tmp[0][0], tmp[-1][1]] = {}
            projclu[c][tmp[0][0], tmp[-1][1]] = {}
            for p in samples:
                projbbc[c][tmp[0][0], tmp[-1][1]][p] = {}
                projclu[c][tmp[0][0], tmp[-1][1]][p] = {}
                projbbc[c][tmp[0][0], tmp[-1][1]][p]['RDR'] = sum(bbc[c][b][p]['RDR'] for b in tmp) / float(len(tmp))
                projbbc[c][tmp[0][0], tmp[-1][1]][p]['BAF'] = sum(bbc[c][b][p]['BAF'] for b in tmp) / float(len(tmp))
            projclu[c][tmp[0][0], tmp[-1][1]] = argmax(dict(Counter([clusters[c][s] for s in tmp])))
            bins = bins[resolution:]
    return projbbc, projclu


def select(bbc, clusters, args):
    alls = set(clusters[c][s] for c in clusters for s in clusters[c])
    count = {idx : {'SIZE' : 0.0, 'CHRS' : set()} for idx in alls}
    totsize = sum(1.0 for c in bbc for s in bbc[c])
    for c in bbc:
        for s in bbc[c]:
            count[clusters[c][s]]['SIZE'] += 1.0
            count[clusters[c][s]]['CHRS'].add(c)

    sel = set(alls)
    if args['st'] is not None:
        sel = set(idx for idx in sel if float(count[idx]['SIZE'] / totsize) >= args['st'])
    if args['ct'] is not None:
        sel = set(idx for idx in sel if len(count[idx]['CHRS']) >= args['ct'])
    s = ['{}:\tSIZE= {},\t# CHRS= {}'.format(idx, count[idx]['SIZE'], count[idx]['CHRS']) for idx in sel]
    sys.stderr.write(info('## Selected clusters: \n{}\n'.format('\n'.join(s))))

    resclu = {}
    resbbc = {}
    for c in clusters:
        if c not in resclu:
            resclu[c] = {}
        assert c not in resbbc
        resbbc[c] = {}
        for s in clusters[c]:
            if clusters[c][s] in sel:
                resclu[c][s] = clusters[c][s]
                resbbc[c][s] = bbc[c][s]

    return {c : resbbc[c] for c in resbbc if len(resbbc[c]) > 0}, {c : resclu[c] for c in resclu if len(resclu[c]) > 0}


def addchr(pos):
    ymin, ymax = plt.ylim()
    corners = []
    prev = 0
    val = pos[0][0]
    for x, s in enumerate(pos):
        if x != 0 and pos[x-1][0] != pos[x][0]:
            plt.plot((x, x), (0, ymax+0.4), '--b', linewidth=0.2)
            corners.append((prev, x, val))
            prev = x
            val = s[0]
    corners.append((prev, x, val))
    ticks = [(int(float(o[1] + o[0] + 1) / 2.0), o[2]) for o in corners]
    plt.xticks([x[0] for x in ticks], [x[1] for x in ticks], rotation=45, ha='center')
    plt.yticks(rotation=0)


def coordinates(args, g=None):
    if g is None:
        if args['xmin'] is not None:
            plt.xlim(xmin=args['xmin'])
        if args['xmax'] is not None:
            plt.xlim(xmax=args['xmax'])
        if args['ymin'] is not None:
            plt.ylim(ymin=args['ymin'])
        if args['ymax'] is not None:
            plt.ylim(ymax=args['ymax'])
    else:
        g.set(xlim=(args['xmin'], args['xmax']), ylim=(args['ymin'], args['ymax']))


def sortchr(x):
    return int(''.join([d for d in x if d.isdigit()]))

def argmax(d):
    return max(d, key=(lambda x : d[x]))

def argmin(d):
    return min(d, key=(lambda x : d[x]))

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def error(msg):
    return "{}{}{}".format("\033[91m\033[1m", msg, "\033[0m")

def warning(msg):
    return "{}{}{}".format("\033[93m\033[1m", msg, "\033[0m")

def log(msg):
    return "{}{}{}".format("\033[95m\033[1m", msg, "\033[0m")

def info(msg):
    return "{}{}{}".format("\033[96m", msg, "\033[0m")

def debug(msg):
    return "{}{}{}".format("\033[92m", msg, "\033[0m")

if __name__ == '__main__':
    main()
