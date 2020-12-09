#!/usr/bin/python3

import sys
import math
import copy
import numpy as np
from scipy.stats import beta

from . import ProgressBar as pb
from . import Supporting as sp
from .ArgParsing import parse_combbo_args


def main(args=None):
    sp.log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_combbo_args(args)
    np.random.seed(seed=args["seed"])
    sp.log(msg="# Reading and checking the bin count files for computing read-depth ratios\n", level="STEP")
    normalbins, tumorbins, chromosomes, normal1, samples1 = readBINs(normalbins=args["normalbins"], tumorbins=args["tumorbins"])
    sp.log(msg="# Reading and checking the allele count files for computing BAF\n", level="STEP")
    tumorbafs, chromosomes2, samples2, normalbafs, normal2 = readBAFs(tumor=args["tumorbafs"], normal=args["normalbafs"])
    if normal1 != normal2 and normal2 is not None: raise ValueError(sp.error("The name of normal sample is different in bin counts and allele counts!"))
    else: normal = normal1

    if samples1 != samples2: raise ValueError(sp.error("The names of tumor samples are different in bin counts and allele counts!"))
    else: samples = samples1

    totalcounts = None
    if args["totalcounts"] is not None:
        sp.log(msg="# Reading and checking the total read count files\n", level="STEP")
        totalcounts = readTotalCounts(filename=args["totalcounts"], samples=samples, normal=normal)

    sp.log(msg="# Combine the bin and allele counts to obtain BAF and RD for each bin\n", level="STEP")
    result = combine(normalbins=normalbins, tumorbins=tumorbins, tumorbafs=tumorbafs, normalbafs=normalbafs, diploidbaf=args["diploidbaf"], totalcounts=totalcounts, chromosomes=chromosomes, samples=samples, normal=normal, mode=args["mode"], bafsd=args["bafsd"], gamma=args["gamma"], draws=args["bootstrap"], verbose=args["verbose"], disable=args["disable"])

    names = list(samples).sort()
    sys.stdout.write("#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\n")
    for key in sorted(result, key=(lambda x : (sp.numericOrder(x[0]), int(x[1]), int(x[2])))):
        for sample in result[key]:
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key[0], key[1], key[2], sample[0], sample[1], sample[2], sample[3], sample[4], sample[5], sample[6]))


def combine(normalbins, tumorbins, tumorbafs, normalbafs, diploidbaf, totalcounts, chromosomes, samples, normal, mode, draws, bafsd, gamma, verbose=False, disable=False):
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
            
            if normalbafs != None:
                normalover = [normalbafs[bi[0], x[1]] for x in tumorover]
            else:
                normalover = None

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
                records = computeBAFs(mode=mode, partition=tpartition, normal=normalover, diploidbaf=diploidbaf, samples=samples, chro=bi[0], draws=draws, bafsd=bafsd, gamma=gamma)
                parsed = {record[0] : record for record in records}

                res[bi] = [(parsed[sample][0], ratios[sample], snps[sample], cov[sample], parsed[sample][1], parsed[sample][2], parsed[sample][3]) for sample in samples]
            else:
                if verbose and normal_count != 0:
                    sp.log(msg='The bin ({}, {}) in chromosomes "{}" has been discarded because there are no covering SNPs in each tumor or normal sample\n'.format(bi[1], bi[2], bi[0]), level="WARN")
                elif verbose and normal_count == 0:
                    sp.log(msg='The bin ({}, {}) in chromosomes "{}" has been discarded because normal read count is zero\n'.format(bi[1], bi[2], bi[0]), level="WARN")
            if not disable: progress_bar.progress(advance=True, msg="Combine bin ({}, {}) in chromosome {}".format(bi[1], bi[2], bi[0]))
    return res


def computeBAFs(mode, partition, normal, diploidbaf, samples, chro, draws, bafsd, gamma):
    if mode == "MOMENTS":
        return moments(partition, samples, normal, diploidbaf, chro, draws, bafsd)
    elif mode == "MIRROR":
        return mirror(partition, samples, normal, diploidbaf)
    elif mode == "BINOMIAL_TEST":
        return testBinomial(partition, samples, normal, diploidbaf, draws, bafsd, gamma)
    else:
        raise ValueError(sp.error("Specified MODE does not exist!"))


def mirror(partition, samples, normal, diploidbaf):
    alphas = {sample : sum(int(min(x[2], x[3])) for x in partition[sample]) for sample in samples}
    betas = {sample : sum(int(max(x[2], x[3])) for x in partition[sample]) for sample in samples}
    mirrorbaf = {sample : (float(alphas[sample]) / float(alphas[sample]+betas[sample])) if (alphas[sample]+betas[sample])>0 else 0.5 for sample in samples}

    if normal != None:
        normalalpha = sum(int(min(x[0], x[1])) for x in normal)
        normalbeta = sum(int(max(x[0], x[1])) for x in normal)
        normalbaf = float(normalalpha) / float(normalalpha + normalbeta)
        normalize = (lambda m, n, d : min(0.5, float(float(m) / float(n)) * 0.5 ) if n > 0 and (0.5 - m) <= d else m)
        mirrorbaf = {sample : normalize(mirrorbaf[sample], normalbaf, diploidbaf) for sample in samples}
        ab = {sample : splitBAF(mirrorbaf[sample], alphas[sample]+betas[sample]) for sample in samples}
        alphas = {sample : ab[sample][0] for sample in samples}
        betas = {sample : ab[sample][1] for sample in samples}
        
    return [(sample, alphas[sample], betas[sample], mirrorbaf[sample]) for sample in samples]


def testBinomial(partition, samples, normal, diploidbaf, draws, bafsd, gamma):
    mirrorres = mirror(partition, samples, normal, diploidbaf)
    return [(record[0], record[1], record[2], 0.5) if isNeutral(record[1], record[2], gamma) else record for record in mirrorres]


def isNeutral(countA, countB, gamma):
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], countA + 1, countB + 1)
    return c_lower <= 0.5 <= c_upper

# def isNeutral(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='beta')
#     return lb <= 0.5 <= ub

# def isNeutral(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='jeffreys')
#     return lb <= 0.5 <= ub


def moments(partition, samples, normal, diploidbaf, chro, draws, bafsd):
    snps = partition
    if normal != None:
        normalsnps = [x for x in normal]
    if draws > 0:
        snps = bootstrap(partition=partition, samples=samples, draws=draws, bafsd=bafsd)
        if normal != None:
            normalsnps = normalbootstrap(points=normal, draws=draws, bafsd=bafsd)

    # Scale BAFs if normal BAFs are given
    Var = (lambda value, mean, size : float(value - mean)**2 / float(size - 1) if size > 1 else 0.0)

    sample_mean = {sample : min(0.5, float(sum(x[4] for x in snps[sample])) / float(len(snps[sample]))) for sample in samples}
    sample_var = {sample : float(sum(Var(x[4], sample_mean[sample], len(snps[sample])) for x in snps[sample])) for sample in samples}

    if normal != None:
        normal_mean = float(sum(x[2] for x in normalsnps )) / float(len(normalsnps))
        norm = float(normal_mean) if normal_mean > 0 else 0.5
        normalize = (lambda m, n, d : min(0.5, (m / n) * 0.5) if (0.5 - m) <= d else m )
        sample_mean = {sample : normalize(sample_mean[sample], norm, diploidbaf) for sample in samples}

    posvar = set(sample for sample in samples if sample_var[sample] > 0.0)
    common = {sample : float( (sample_mean[sample] * (1.0 - sample_mean[sample])) / sample_var[sample]) - 1.0 for sample in posvar}
    alphas = {sample : float(sample_mean[sample] * common[sample]) if sample in posvar else sum(min(x[2], x[3]) for x in snps[sample]) for sample in samples}
    betas = {sample : float((1.0 - sample_mean[sample]) * common[sample]) if sample in posvar else sum(max(x[2], x[3]) for x in snps[sample]) for sample in samples}
    roundings = {sample : roundAlphasBetas(sample_mean[sample], alphas[sample], betas[sample]) for sample in samples}
    
    return [(sample, roundings[sample][0], roundings[sample][1], sample_mean[sample]) for sample in samples]


def bootstrap(partition, samples, draws, bafsd):
    boot = {sample : [] for sample in samples}
    gauss = np.random.normal
    binomial = np.random.binomial
    poisson = np.random.poisson
    avg_cov = {sample : float(sum(x[2]+x[3] for x in partition[sample])) / float(len(partition[sample])) for sample in samples}
    for sample in samples:
        for snp in partition[sample]:
            # OPTION 2
            # bafs = [snp[4] for i in range(draws)]
            # covs = map(float, list(poisson(lam=avg_cov[snp[0]], size=draws)))
            # covs = [c if c > 0.0 else float(avg_cov[snp[0]]) for c in covs]
            # mid = int(len(covs)/2)
            # alleles = [binomial(n=cov, p=baf) for cov, baf in zip(covs[:mid], bafs[:mid])]
            # alleles += [binomial(n=cov, p=(1.0 - baf)) for cov, baf in zip(covs[mid:], bafs[mid:])]
            # boot[sample].extend([(snp[0], snp[1], int(alleles[i]), int(covs[i]-alleles[i]), float(min(alleles[i], covs[i]-alleles[i])) / float(covs[i]) ) for i in range(draws)])
            covs = map(float, list(poisson(lam=avg_cov[snp[0]], size=draws)))
            covs = [c if c > 0.0 else float(avg_cov[snp[0]]) for c in covs]
            ebafs = sorted([gauss(snp[4], bafsd) for cov in covs])
            winsorized = int(len(ebafs) * 0.1)
            ebafs = [ebafs[winsorized] for i in range(winsorized)] + ebafs[winsorized:-winsorized] + [ebafs[-winsorized-1] for i in range(winsorized)]
            assert(len(ebafs) == len(covs))
            alleles = [splitBAF(ebaf, cov) for (ebaf, cov) in zip(ebafs, covs)]
            boot[sample].extend([(snp[0], snp[1], alleles[i][0], alleles[i][1], float(min(alleles[i][0], alleles[i][1])) / float(alleles[i][0] + alleles[i][1]) ) for i in range(draws)])
    return {sample : partition[sample] + boot[sample] for sample in samples}


def normalbootstrap(points, draws, bafsd):
    boot = []
    gauss = np.random.normal
    binomial = np.random.binomial
    poisson = np.random.poisson
    avg_cov = float(sum(int(x[0]) + int(x[1]) for x in points)) / float(len(points))
    for snp in points:
        # OPTION 2:
        # bafs = [snp[2] for i in range(draws)]
        # covs = map(float, list(poisson(lam=avg_cov, size=draws)))
        # covs = [c if c > 0.0 else float(avg_cov) for c in covs]
        # mid = int(len(covs)/2)
        # alleles = [binomial(n=cov, p=baf) for cov, baf in zip(covs[:mid], bafs[:mid])]
        # alleles += [binomial(n=cov, p=(1.0 - baf)) for cov, baf in zip(covs[mid:], bafs[mid:])]
        # boot.extend([(int(alleles[i]), int(covs[i]-alleles[i]), float(min(alleles[i], covs[i]-alleles[i])) / float(covs[i]) ) for i in range(draws)])
        covs = map(float, list(poisson(lam=avg_cov, size=draws)))
        covs = [c if c > 0.0 else float(avg_cov) for c in covs]
        ebafs = [gauss(snp[2], bafsd) for cov in covs]
        winsorized = int(len(ebafs) * 0.1)
        ebafs = [ebafs[winsorized] for i in range(winsorized)] + ebafs[winsorized:-winsorized] + [ebafs[-winsorized-1] for i in range(winsorized)]
        assert(len(ebafs) == len(covs))
        alleles = [splitBAF(ebaf, cov) for (ebaf, cov) in zip(ebafs, covs)]
        boot.extend([(int(alleles[i][0]), int(alleles[i][1]), float(min(alleles[i][0], alleles[i][1])) / float(alleles[i][0]+alleles[i][1]) ) for i in range(draws)])
    return points + boot
    

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


def readBAFs(tumor, normal=None):
    tumorBAFs = {}
    normal_chr = set()
    tumor_chr = set()

    if normal != None:
        normalBAFs = {}
        normal_sample = set()
        normal_chr = set()

        with open(normal, 'r') as f:
            for line in f:
                parsed = line.strip().split()[:5]
                normal_chr.add(parsed[1])
                normal_sample.add(parsed[0])
                pos = int(parsed[2])
                ref = int(parsed[3])
                alt = int(parsed[4])
                if (parsed[1], pos) not in normalBAFs:
                    normalBAFs[parsed[1], pos] = (ref, alt, float(min(ref, alt)) / float(ref + alt))
                else:
                    raise ValueError(sp.error("A position is present multiple times in the normal sample!"))

        if len(normal_sample) != 1:
            raise ValueError("Found multiple samples in the normal BAF file!")
        if not tumor_chr <= normal_chr:
            raise ValueError("The chromosomes in tumor are not present in the normal!")

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
            if normal != None and (parsed[1], pos) not in normalBAFs:
                raise ValueError("The SNP at position {} in chromosome {} is covered in a tumor sample but not in the normal!".format(pos, parsed[1]))
            # if normal is not None:
            #     try:
            #         scale = ref + alt
            #         normalbaf = normalBAFs[parsed[1], pos][2]
            #         baf = float(baf) / float(normalbaf) * float(0.5) if normalbaf > 0 else 0.5
            #         #baf = min(baf, 1.0 - baf)
            #         baf = min(baf, 0.5)
            #         alt, ref = splitBAF(baf, scale)
            #     except KeyError:
            #         raise ValueError("The SNP at position {} in chromosome {} is covered in a tumor sample but not in the normal!".format(pos, parsed[1]))
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

    if normal is None:
        return tumorBAFs, chromosomes, samples, None, None
    else:
        return tumorBAFs, chromosomes, samples, normalBAFs, normal_sample.pop()


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


def roundAlphasBetas(baf, alpha, beta):
    BAF = float(baf)
    BAF = min(BAF, 1.0 - BAF)
    ALPHA = min(alpha, beta)
    BETA = max(alpha, beta)

    roundings = []
    roundings.append((int(math.floor(ALPHA)), int(math.floor(BETA))))
    roundings.append((int(math.floor(ALPHA)), int(math.ceil(BETA))))
    roundings.append((int(math.ceil(ALPHA)), int(math.floor(BETA))))
    roundings.append((int(math.ceil(ALPHA)), int(math.ceil(BETA))))
    roundings = [(int(min(a,b)), int(max(a,b))) for (a, b) in roundings]

    estimations = [float(a) / float(a+b) if a+b>0 else 1.0 for (a, b) in roundings]
    diff = [abs(est - BAF) for est in estimations]
    return roundings[np.argmin(diff)]


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
