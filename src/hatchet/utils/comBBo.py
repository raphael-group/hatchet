#!/usr/bin/python2

import sys
import math
import copy
import numpy as np
from scipy.stats import beta

import ProgressBar as pb
import Supporting as sp
from ArgParsing import parse_combbo_args

from collections import defaultdict
from collections import deque



def main(args=None):
    sp.log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_combbo_args(args)
    sp.logArgs(args, 80)
    np.random.seed(seed=args["seed"])
    sp.log(msg="# Reading and checking the bin count files for computing read-depth ratios\n", level="STEP")
    normalbins, tumorbins, chromosomes, normal, samples1 = readBINs(normalbins=args["normalbins"], tumorbins=args["tumorbins"])
    sp.log(msg="# Reading and checking the allele count files for computing BAF\n", level="STEP")
    tumorbafs, chromosomes2, samples2 = readBAFs(tumor=args["tumorbafs"])
    if samples1 != samples2: raise ValueError(sp.error("The names of tumor samples are different in bin counts and allele counts!"))
    else: samples = samples1
    if args['phase'] is not None:
        sp.log(msg="# Reading phases of heterozygous germline SNPs\n", level="STEP")
        phase = readPhase(args["phase"])
    else:
        phase = None

    totalcounts = None
    if args["totalcounts"] is not None:
        sp.log(msg="# Reading and checking the total read count files\n", level="STEP")
        totalcounts = readTotalCounts(filename=args["totalcounts"], samples=samples, normal=normal)

    sp.log(msg="# Combine the bin and allele counts to obtain BAF and RD for each bin\n", level="STEP")
    result = combine(normalbins=normalbins, tumorbins=tumorbins, tumorbafs=tumorbafs, diploidbaf=args["diploidbaf"], 
                     totalcounts=totalcounts, chromosomes=chromosomes, samples=samples, normal=normal, gamma=args["gamma"], 
                     verbose=args["verbose"], disable=args["disable"], phase=phase, block=args["block"])

    names = list(samples).sort()
    sys.stdout.write("#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\n")
    for key in sorted(result, key=(lambda x : (sp.numericOrder(x[0]), int(x[1]), int(x[2])))):
        for sample in result[key]:
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key[0], key[1], key[2], sample[0], sample[1], sample[2], sample[3], sample[4], sample[5], sample[6]))


def combine(normalbins, tumorbins, tumorbafs, diploidbaf, totalcounts, chromosomes, samples, normal, gamma, verbose=False, disable=False, phase=None, block=None):
    res = {}
    ctumorbafs = {c : sorted(tumorbafs[c], key=(lambda x : x[1])) for c in tumorbafs}
    if not disable: progress_bar = pb.ProgressBar(total=len(tumorbins), length=40, verbose=verbose)
    for bi in sorted(tumorbins, key=(lambda x : (sp.numericOrder(x[0]), int(x[1]), int(x[2])))):
        if bi[0] in ctumorbafs:
            # Extract counts
            normal_count = normalbins[bi][1]
            tumor_counts = {x[0] : x[1] for x in tumorbins[bi]}

            # Extract overlapping SNPs: option 2
            #tumorover = [x for x in ctumorbafs[bi[0]] if bi[1] <= x[1] <= bi[2]]
            #ctumorbafs[bi[0]] = [x for x in ctumorbafs[bi[0]] if x not in tumorover]

            tumorover = []
            for i, x in enumerate(ctumorbafs[bi[0]]):
                if bi[1] <= x[1] <= bi[2]:
                    tumorover.append(x)
                elif x[1] > bi[2]:
                    ctumorbafs[bi[0]] = ctumorbafs[bi[0]][i:]
                    break
            
            # Partition the overlapping SNPs by samples
            tpartition = {sample : [x for x in tumorover if x[0] == sample] for sample in samples}

            # Check the bin to be covered in each sample and non-zero count in normal
            if sum(len(tpartition[key]) != 0 for key in tpartition) == len(samples) and normal_count != 0:
                # Compute ratios
                ratios = {sample : float(tumor_counts[sample]) / float(normal_count) for sample in samples}

                # Compute normalizing factor by total number of reads if provided
                if totalcounts is not None:
                    totalfactor = {sample : float(totalcounts[normal]) / float(totalcounts[sample])  for sample in samples}
                    ratios = {sample : float(ratios[sample]) * float(totalfactor[sample]) for sample in samples}

                # Compute number of SNPs covering each bin and the average coverage
                snps = {sample : len(tpartition[sample]) for sample in samples}
                cov = {sample : float(sum(x[2]+x[3] for x in tpartition[sample])) / float(len(tpartition[sample])) for sample in samples}

                #Compute BAFs
                if phase is None:
                    records = computeBAFs(partition=tpartition, diploidbaf=diploidbaf, samples=samples)
                else:
                    records = computeBAFs(partition=tpartition, diploidbaf=diploidbaf, samples=samples, phase=phase[bi[0]], block=block)
                parsed = {record[0] : record for record in records}

                res[bi] = [(parsed[sample][0], ratios[sample], snps[sample], cov[sample], parsed[sample][1], parsed[sample][2], parsed[sample][3]) for sample in samples]
            else:
                if verbose and normal_count != 0:
                    sp.log(msg='The bin ({}, {}) in chromosomes "{}" has been discarded because there are no covering SNPs in each tumor or normal sample\n'.format(bi[1], bi[2], bi[0]), level="WARN")
                elif verbose and normal_count == 0:
                    sp.log(msg='The bin ({}, {}) in chromosomes "{}" has been discarded because normal read count is zero\n'.format(bi[1], bi[2], bi[0]), level="WARN")
        if not disable: progress_bar.progress(advance=True, msg="Combine bin ({}, {}) in chromosome {}".format(bi[1], bi[2], bi[0]))
    return res


def computeBAFs(partition, samples, diploidbaf, phase=None, block=0):
    if phase is None:
        tpartition = partition
    else:
        select = (lambda L, s : blocking(filter(lambda o : o[1] in phase, L), s, phase, block))
        tpartition = {sample : select(partition[sample], sample) for sample in samples}

    alphas = {sample : sum(int(min(x[2], x[3])) for x in tpartition[sample]) for sample in samples}
    betas = {sample : sum(int(max(x[2], x[3])) for x in tpartition[sample]) for sample in samples}
    mirrorbaf = {sample : (float(alphas[sample]) / float(alphas[sample]+betas[sample])) if (alphas[sample]+betas[sample])>0 else 0.5 for sample in samples}
        
    return [(sample, alphas[sample], betas[sample], mirrorbaf[sample]) for sample in samples]


def blocking(L, sample, phase, blocksize):
    result = []
    if len(L) == 0:
        return result
    que = deque(sorted(L, key=(lambda v : v[1])))
    omap = {}
    blocks = {}
    for bk in range(min(o[1] for o in L), max(o[1] for o in L) + 1, blocksize):
        block = (sample, bk, 0, 0)
        while que and bk <= que[0][1] < bk + blocksize:
            o = que.popleft()
            if phase[o[1]] == '0|1':
                block = (sample, bk, block[2] + o[2], block[3] + o[3])
            elif phase[o[1]] == '1|0':
                block = (sample, bk, block[2] + o[3], block[3] + o[2])
            else:
                assert False, 'Found a wrong phase value'
            omap[o] = bk
        if block[2] + block[3] > 0:
            result.append(block)
    return result


def readBINs(normalbins, tumorbins):
    normalBINs = {}
    tumorBINs = {}
    normal = set()
    samples = set()
    normal_chr = set()
    tumor_chr = set()

    # Read normal bin counts
    with open(normalbins, 'r') as f:
        for line in f:
            parsed = line.strip().split()[:5]
            normal_chr.add(parsed[1])
            normal.add(parsed[0])
            if (parsed[1], int(parsed[2]), int(parsed[3])) not in normalBINs:
                normalBINs[parsed[1], int(parsed[2]), int(parsed[3])] = (parsed[0], int(parsed[4]))
            else:
                raise ValueError(sp.error("Found multiple lines for the same interval in the normal bin counts!"))

    # Check normal bin counts
    if len(normal) > 1:
        raise ValueError(sp.error("Found multiple samples in normal bin counts!"))
    prev_r = -1
    prev_c = -1
    for key in sorted(normalBINs, key=(lambda x : (sp.numericOrder(x[0]), int(x[1]), int(x[2])))):
        l, r = int(key[1]), int(key[2])
        if l > r and prev_c == key[0]:
            raise ValueError(sp.error("Found an interval with START {} greater than END {} in normal bin counts!".format(key[1], key[2])))
        if l < prev_r and prev_c == key[0]:
            raise ValueError(sp.error("Found overlapping intervals one ending with {} and the next starting with {} in normal bin counts!".format(prev_r, key[1])))
        prev_r = r
        prev_c = key[0]

    # Read tumor bin counts
    with open(tumorbins, 'r') as f:
        for line in f:
            parsed = line.strip().split()[:5]
            tumor_chr.add(parsed[1])
            samples.add(parsed[0])
            try:
                tumorBINs[parsed[1], int(parsed[2]), int(parsed[3])].add((parsed[0], int(parsed[4])))
            except KeyError:
                tumorBINs[parsed[1], int(parsed[2]), int(parsed[3])] = set()
                tumorBINs[parsed[1], int(parsed[2]), int(parsed[3])].add((parsed[0], int(parsed[4])))

    # Check tumor bin counts
    prev_r = -1
    prev_c = -1
    num_samples = len(samples)
    for key in sorted(tumorBINs, key=(lambda x : (sp.numericOrder(x[0]), int(x[1]), int(x[2])))):
        l, r = int(key[1]), int(key[2])
        if len(tumorBINs[key]) != num_samples:
            raise ValueError(sp.error("Found multiple lines for the same interval in the tumor bin counts!"))
        if l > r and prev_c == key[0]:
            raise ValueError(sp.error("Found an interval with START {} greater than END {} in tumor bin counts!".format(key[1], key[2])))
        if l < prev_r and prev_c == key[0]:
            raise ValueError(sp.error("Found overlapping intervals one ending with {} and the next starting with {} in tumor bin counts!".format(prev_r, key[1])))
        prev_r = r
        prev_c = key[0]

    if normal_chr != tumor_chr:
        raise ValueError(sp.error("The chromosomes in normal and tumor bin counts are different!"))
    if set(normalBINs) != set(tumorBINs):
        raise ValueError(sp.error("The bins of the normal and tumor samples are different!"))

    chromosomes = sorted(list(normal_chr), key=sp.numericOrder)

    return normalBINs, tumorBINs, chromosomes, normal.pop(), samples


def readBAFs(tumor):
    tumorBAFs = {}
    tumor_chr = set()

    # Read tumor bafs
    samples = set()
    with open(tumor, 'r') as f:
        for line in f:
            parsed = line.strip().split()[:5]
            sample = parsed[0]
            chromosome = parsed[1]
            pos = int(parsed[2])
            ref = int(parsed[3])
            alt = int(parsed[4])
            tumor_chr.add(chromosome)
            samples.add(sample)
            baf = float(min(ref, alt)) / float(ref+alt) if ref+alt > 0 else 0.5
            try:
                tumorBAFs[parsed[1]].append((parsed[0], pos, ref, alt, baf))
            except KeyError:
                tumorBAFs[parsed[1]] = [(parsed[0], pos, ref, alt, baf)]

    # Check tumor bafs
    for key in tumorBAFs:
        tumorBAFs[key].sort(key=(lambda x : x[1]))
        if len(tumorBAFs[key]) > len(set((x[0], x[1]) for x in tumorBAFs[key])):
            raise ValueError(sp.error("A position is present multiple times in the tumor samples!"))

    chromosomes = sorted(list(tumor_chr), key=sp.numericOrder)

    return tumorBAFs, chromosomes, samples


def readPhase(f):
    phased = defaultdict(lambda : dict())
    with open(f, 'r') as i:
        for l in i:
            p = l.strip().split()
            if len(l) > 1 and p[0][0] != '#':
                zeroone = '0|1' in l
                onezero = '1|0' in l
                if zeroone or onezero:
                    if zeroone and onezero:
                        raise ValueError('Found a record in phased positions which contains both phases 0|1 and 1|0!')
                    if p[0] in phased[p[0]]:
                        raise ValueError('Found a duplicate phased position!')
                    phased[p[0]][int(p[1])] = '0|1' if zeroone else '1|0'
    return phased


def splitBAF(baf, scale):
    BAF = float(baf)
    BAF = min(BAF, 1.0 - BAF)
    SUM = float(scale)

    roundings = []
    roundings.append((int(math.floor(BAF * SUM)), int(math.floor((1.0 - BAF) * SUM))))
    roundings.append((int(math.floor(BAF * SUM)), int(math.ceil((1.0 - BAF) * SUM))))
    roundings.append((int(math.ceil(BAF * SUM)), int(math.floor((1.0 - BAF) * SUM))))
    roundings.append((int(math.ceil(BAF * SUM)), int(math.ceil((1.0 - BAF) * SUM))))
    roundings = [(int(min(a,b)), int(max(a,b))) for (a, b) in roundings]

    estimations = [float(a) / float(a+b) if a+b>0 else 1.0 for (a, b) in roundings]
    diff = [abs(est - BAF) for est in estimations]
    best = np.argmin(diff)
    return roundings[best][0], roundings[best][1]


def readTotalCounts(filename, samples, normal):
    normalfound = False
    counts = {}
    found = set()
    with open(filename, 'r') as f:
        for line in f:
            parsed = line.strip().split()
            if parsed[0] in found:
                raise ValueError(sp.error("Found multiple total read counts for the same sample {}!".format(parsed[0])))
            if parsed[0] == normal:
                normalfound = True
            else:
                found.add(parsed[0])
            counts[parsed[0]] = int(parsed[1])
    if samples < found:
        raise ValueError(sp.error("Found total read counts for samples that are not present in the input!"))
    elif found < samples:
        raise ValueError(sp.error("Missing total read counts for some samples present in the input!"))
    elif not normalfound:
        raise ValueError(sp.error("Missing total read counts for normal sample!"))
    else:
        return counts


if __name__ == '__main__':
    main()
