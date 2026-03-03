import sys
import os
import argparse
import math
import numpy as np

from collections import deque
from hatchet import __version__


def parse_args():
    description = (
        "This method takes in input multiple samples from the same patient, where each sample is a "
        "segmented CNV file produced by GATK4 CNV pipeline, and produces a BB input file for HATCHet2."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "INPUT",
        type=str,
        help=(
            "A white-space-separated list between apices where each element is a segmented CNV file produced by "
            "GATK4 CNV pipeline. The file format is describe in the HATCHet2's repository."
        ),
    )
    parser.add_argument(
        "--samples",
        type=str,
        required=False,
        default=None,
        help=(
            "A white-space-separated list containing the name of sample in the same order as given (default: "
            "filenames are used)."
        ),
    )
    parser.add_argument(
        "-b",
        "--binsize",
        required=False,
        default="50kb",
        type=str,
        help='Size of the bins, specified as a full number or using the notations either "kb" or "Mb" (default: 50kb).',
    )
    parser.add_argument(
        "-r",
        "--devRDR",
        type=float,
        required=False,
        default=0.05,
        help="Standard deviation for the RDR of bins obtained from RDR of the corresponding segment (default: 0,05).",
    )
    parser.add_argument(
        "-a",
        "--devBAF",
        type=float,
        required=False,
        default=0.02,
        help="Standard deviation for the RDR of bins obtained from RDR of the corresponding segment (default: 0.02).",
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        required=False,
        default=None,
        help="Starting seed for random number generator (default: not specified).",
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()

    samples = args.INPUT.strip().split()
    for s in samples:
        if not os.path.isfile(s):
            raise ValueError("ERROR: sample {} does not exist!".format(s))

    if args.samples:
        names = args.samples.strip().split()
        if len(names) != len(samples):
            raise ValueError(
                "ERROR: a sample name must be specified for every given sample!"
            )
        samples = {n: s for n, s in zip(names, samples)}
    else:
        samples = {n: s for n, s in zip(samples, samples)}

    size = 0
    try:
        if args.binsize[-2:] == "kb":
            size = int(args.binsize[:-2]) * 1000
        elif args.binsize[-2:] == "Mb":
            size = int(args.binsize[:-2]) * 1000000
        else:
            size = int(args.binsize)
    except (IndexError, ValueError):
        raise ValueError(
            'Size must be a number, optionally ending with either "kb" or "Mb"!'
        )

    if args.devRDR < 0.0:
        raise ValueError("ERROR: deviation of RDR must be grater or equal than 0.0!")

    if args.devBAF < 0.0:
        raise ValueError("ERROR: deviation of BAF must be grater or equal than 0.0!")

    if args.seed:
        np.random.seed(args.seed)

    return {
        "samples": samples,
        "binsize": size,
        "devrdr": args.devRDR,
        "devbaf": args.devBAF,
    }


def main():
    log("Parsing and checking input arguments")
    args = parse_args()

    log("Reading input files")
    segs = read_segs(args["samples"])

    log("Segmenting all samples jointly")
    segs = jsegmentation(segs)

    log("Binning the segments across all samples")
    bins = binning(segs, size=args["binsize"], drdr=args["devrdr"], dbaf=args["devbaf"])

    log("Writing the corresponding BB file")
    print(
        "\t".join(
            [
                "#CHR",
                "START",
                "END",
                "SAMPLE",
                "RD",
                "#SNPS",
                "COV",
                "ALPHA",
                "BETA",
                "BAF",
            ]
        )
    )
    orderchrs = lambda x: int("".join([l for l in x if l.isdigit()]))
    nsnp = lambda L: int(round(L / 1000.0))
    cov = lambda R: R * 60.0
    AB = lambda R, B, L: list(splitBAF(baf=B, scale=max(1000, R * L / 250.0)))
    dat = lambda c, b, p, L, R, B: [R, nsnp(L), cov(R)] + AB(R, B, L) + [B]
    row = lambda c, b, p: [c, b[0], b[1], p] + dat(
        c, b, p, float(b[1] - b[0]), bins[c][b][p]["RDR"], bins[c][b][p]["BAF"]
    )
    print(
        "\n".join(
            [
                "\t".join(map(str, row(c, b, p)))
                for c in sorted(bins, key=orderchrs)
                for b in sorted(bins[c])
                for p in sorted(bins[c][b])
            ]
        )
    )


def read_segs(samples):
    def read_sample(s):
        d = {}
        with open(s, "r") as i:
            for l in i:
                if l[0] != "#" and l[0] != "@" and len(l) > 1:
                    p = l.strip().split()
                    if p[0] != "CONTIG":
                        h = p[0]
                        seg = (int(p[1]), int(p[2]))
                        if p[6].upper() not in ["NAN", "NONE"] and p[9].upper() not in [
                            "NAN",
                            "NONE",
                        ]:
                            rdr = math.pow(2.0, float(p[6]))
                            if rdr < 0.0:
                                raise ValueError(
                                    "ERROR: the log of RDR expected to be positive but value {} found!".format(
                                        rdr
                                    )
                                )
                            baf = float(p[9])
                            if not 0.0 <= baf <= 0.5:
                                raise ValueError(
                                    "ERROR: BAF expected to be in [0, 0.5] but value {} found!".format(
                                        baf
                                    )
                                )
                            if h not in d:
                                d[h] = {}
                            assert (
                                seg not in d[h]
                            ), "ERROR: Found a duplicate segment {}:{}-{}".format(
                                h, seg[0], seg[1]
                            )
                            check = lambda b: (b[0] <= b[1] <= seg[0]) or (
                                seg[1] <= b[0] <= b[1]
                            )
                            assert False not in set(
                                check(b) for b in d[h]
                            ), "ERROR: found overlapping segment {}:{}-{}".format(
                                h, seg[0], seg[1]
                            )
                            d[h][seg] = {"RDR": rdr, "BAF": baf}
        return d

    return {s: read_sample(samples[s]) for s in samples}


def jsegmentation(segs):
    chrs = set(c for p in segs for c in segs[p])
    bk = lambda p, c: set(k for s in segs[p][c] for k in s)
    bks = {
        c: sorted(set(k for p in segs if c in segs[p] for k in bk(p, c))) for c in chrs
    }
    counts = {c: {s: 0 for s in zip(bks[c][:-1], bks[c][1:])} for c in bks}
    bmap = {p: {c: {s: None for s in counts[c]} for c in counts} for p in segs}

    for p in segs:
        for c in segs[p]:
            bk = deque(bks[c])
            left = -1
            right = bk.popleft()
            for l, r in sorted(segs[p][c], key=(lambda x: x[0])):
                while right != r:
                    left = right
                    right = bk.popleft()
                    if l <= left and right <= r:
                        counts[c][left, right] += 1
                        bmap[p][c][left, right] = (l, r)

    tak = {c: set(s for s in counts[c] if counts[c][s] == len(segs)) for c in counts}
    tak = {c: tak[c] for c in tak if len(tak[c]) > 0}

    tot = float(sum(s[1] - s[0] for c in counts for s in counts[c]))
    cov = sum(s[1] - s[0] for c in tak for s in tak[c])
    log(
        "##Coverage from joint segmentation is {0:.2f}%".format(
            float(100.0 * cov) / tot
        )
    )

    res = {
        p: {c: {s: segs[p][c][bmap[p][c][s]] for s in tak[c]} for c in tak}
        for p in segs
    }

    chrs = set(c for p in res for c in res[p])
    pos = {c: set(s for p in res for s in res[p][c]) for c in chrs}

    join = {p: {c: {} for c in segs[p]} for p in segs}
    for c in chrs:
        l = None
        r = None
        val = None
        for s in sorted(pos[c], key=(lambda x: x[0])):
            if l is None:
                l = s[0]
                val = {p: res[p][c][s] for p in res}
            elif val != {p: res[p][c][s] for p in res}:
                for p in join:
                    join[p][c][l, r] = val[p]
                l = s[0]
                val = {p: res[p][c][s] for p in res}
            r = s[1]
        for p in join:
            join[p][c][l, r] = val[p]

    assert False not in set(
        cov == sum(s[1] - s[0] for c in join[p] for s in join[p][c]) for p in res
    )

    return join


def binning(segs, size, drdr, dbaf):
    tg = list(segs.keys())[0]
    chrs = segs[tg].keys()
    pos = {c: sorted([s for s in segs[tg][c]], key=(lambda x: x[0])) for c in chrs}
    norm = np.random.normal
    BAF = (
        lambda x: min(x, 1.0 - x)
        if 0.0 <= x <= 1.0
        else (BAF(1.0) if x > 1.0 else BAF(-x))
    )
    gen = lambda d: {
        "RDR": norm(d["RDR"], drdr),
        "BAF": BAF(norm(d["BAF"], dbaf)),
    }

    bins = {}
    for c in chrs:
        bins[c] = {}
        for s in pos[c]:
            assert (
                s[1] - s[0] > 0
            ), "ERROR: START and END cannot be equal: {}:{}-{}".format(c, s[0], s[1])
            part = list(range(s[0], s[1], size))
            part = part + [s[1]] if part[-1] != s[1] else part
            for b in zip(part[:-1], part[1:]):
                bins[c][b] = {p: gen(segs[p][c][s]) for p in segs}

    return bins


def splitBAF(baf, scale):
    BAF = float(baf)
    assert 0.0 <= BAF <= 0.5
    SUM = float(scale)

    roundings = []
    roundings.append((int(math.floor(BAF * SUM)), int(math.floor((1.0 - BAF) * SUM))))
    roundings.append((int(math.floor(BAF * SUM)), int(math.ceil((1.0 - BAF) * SUM))))
    roundings.append((int(math.ceil(BAF * SUM)), int(math.floor((1.0 - BAF) * SUM))))
    roundings.append((int(math.ceil(BAF * SUM)), int(math.ceil((1.0 - BAF) * SUM))))
    roundings = [(int(min(a, b)), int(max(a, b))) for (a, b) in roundings]

    estimations = [
        float(a) / float(a + b) if a + b > 0 else 1.0 for (a, b) in roundings
    ]
    diff = [abs(est - BAF) for est in estimations]
    best = np.argmin(diff)

    return roundings[best][0], roundings[best][1]


def log(M):
    sys.stderr.write("#" + M + "\n")


if __name__ == "__main__":
    main()
