#!/usr/bin/python2

import os, shutil
import sys
import math
import copy
import numpy as np

from ArgParsing import parse_clubb_args
import Supporting as sp


def main():
    sp.log(msg="# Parsing and checking input arguments\n", level="STEP")
    args = parse_clubb_args()

    sp.log(msg="# Reading the combined BB file\n", level="STEP")
    combo, samples = readBB(args["bbfile"])

    sp.log(msg="# Format data to cluster\n", level="STEP")
    points, bintoidx = getPoints(data=combo, samples=samples)

    clouds = None
    if args["cloud"] > 0 :
        sp.log(msg="# Bootstrap each bin for clustering\n", level="STEP")
        clouds = generateClouds(points=points, density=args["cloud"], seed=args["seed"], sdeven=args["ratiodeviation"], sdodd=args["bafdeviation"])

    sp.log(msg="# Clustering bins by RD and BAF across tumor samples\n", level="STEP")
    mus, sigmas, clusterAssignments, numPoints, numClusters = cluster(points=points, output=args["outsegments"], samples=samples, clouds=clouds, K=args["initclusters"], sf=args["tuning"], bnpydir=args["bnpydir"])

    if args['rdtol'] > 0.0 or args['baftol'] > 0.0:
        sp.log(msg="# Refining clustering using given tolerances\n", level="STEP")
        before = len(set(clusterAssignments))
        clusterAssignments, numClusters = refineClustering(combo=combo, assign=clusterAssignments, assignidx=bintoidx, samples=samples, rdtol=args['rdtol'], baftol=args['baftol'])
        sp.log(msg='The number of clusters have been reduced from {} to {} with given tolerances\n'.format(before, numClusters), level='INFO')

    names = list(samples).sort()

    sp.log(msg="# Writing BBC output with resulting clusters\n", level="STEP")
    if args["outbins"] is None: outbins = sys.stdout
    else: outbins = open(args["outbins"], 'w')
    outbins.write("#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\tCLUSTER\n")
    for key in sorted(combo, key=(lambda x : (sp.numericOrder(x[0]), int(x[1]), int(x[2])))):
        for sample in combo[key]:
            outbins.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key[0], key[1], key[2], sample[0], sample[1], sample[2], sample[3], sample[4], sample[5], sample[6], clusterAssignments[bintoidx[key]]))

    sp.log(msg="# Segmenting bins\n", level="STEP")
    clusters = {cluster : set(key for key in combo if clusterAssignments[bintoidx[key]] == cluster) for cluster in set(clusterAssignments)}
    segments = segmentBins(bb=combo, clusters=clusters, samples=samples)

    if args["diploidbaf"] != None:
        sp.log(msg="# Determining the largest cluster as diploid or tetraploid and rescaling all the clusters inside the threshold accordingly\n", level="STEP")
        segments = scaleBAF(segments=segments, samples=samples, diploidbaf=args["diploidbaf"])

    sp.log(msg="# Writing REF output with resulting segments\n", level="STEP")
    if args["outsegments"] is None: outsegments = sys.stdout
    else: outsegments = open(args["outsegments"], 'w')
    outsegments.write("#ID\tSAMPLE\t#BINS\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\n")
    for key in sorted(segments):
        for sample in segments[key]:
            record = segments[key][sample]
            outsegments.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, sample, record[0], record[1], record[2], record[3], record[4], record[5], record[6]))


def readBB(bbfile):
    read = {}
    samples = set()
    with open(bbfile, 'r') as f:
        for line in f:
            if line[0] != '#':
                parsed = line.strip().split()
                chromosome = parsed[0]
                start = int(parsed[1])
                end = int(parsed[2])
                sample = parsed[3]
                rd = float(parsed[4])
                snps = int(parsed[5])
                cov = float(parsed[6])
                alpha = int(parsed[7])
                beta = int(parsed[8])
                baf = float(parsed[9])

                samples.add(sample)
                try:
                    read[chromosome, start, end].append((sample, rd, snps, cov, alpha, beta, baf))
                except KeyError:
                    read[chromosome, start, end] = [(sample, rd, snps, cov, alpha, beta, baf)]
    return read, samples


def getPoints(data, samples):
    idx = 0
    points = []
    bintoidx = {}
    for bi in sorted(data, key=(lambda x : (sp.numericOrder(x[0]), int(x[1]), int(x[2])))):
        partition = {}
        for x in data[bi]:
            if x[0] in partition:
                raise ValueError(sp.error("Found a bin ({}, {}) in chromosome {} defined multiple times for the same sample!".format(bi[1]. bi[2], bi[0])))
            else:
                partition[x[0]] = [x[1], x[-1]]
        if len(partition) != len(samples):
            raise ValueError(sp.error("Found a bin ({}, {}) in chromosome {} that is not covered in all the samples!".format(bi[1]. bi[2], bi[0])))
        points.append([e for sample in samples for e in partition[sample]])
        bintoidx[bi] = idx
        idx += 1
    return points, bintoidx


def cluster(points, output, samples, clouds=None, K=15, sf=0.01, bnpydir=None):
    """
    Clusters a set of data points lying in an arbitrary number of clusters.
    Arguments:
        data (list of lists of floats): list of data points to be clustered.
        sampleName (string): The name of the input sample.
        sf (float): Tuning parameter for clustering; used to determine initial size of
                    distribution covariances. Small sf indicates a belief that clusters
                    are of small size.
    Returns:
        mus (list of lists of floats): List of cluster means.
        sigmas (list of 2D lists of floats): List of cluster covariances.
        clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
                                            j at index i means the ith interval has been assigned to the
                                            jth meta-interval.
        numPoints (list of ints): Number of points assigned to each cluster
        numClusters (int): The number of clusters.
    """
    sp.log(msg="## Loading BNPY\n", level="INFO")
    tmp = os.path.splitext(output)[0] + "_TMPDIR/"
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    os.environ["BNPYOUTDIR"] = tmp
    sys.path.append(bnpydir)
    import bnpy

    sp.log(msg="## Clustering...\n", level="INFO")
    total = list(points)
    if clouds is not None:
        total.extend(list(clouds))
    npArray = np.array(total)
    Data = bnpy.data.XData(X=npArray)
    Data.name = "Clustering tumor samples by RD and BAF"
    Data.summary = "Clustering the following samples: {}".format(",".join(samples))

    ##K = 15
    if Data.X.shape[0] < K:
	    K = Data.X.shape[0]

#        hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=100, nTask=1, K=K, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=sf, doWriteStdOut=False)
#        hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', moves='birth,merge', sF=0.0001, doWriteStdOut=False)
#        hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nTask=1, K=K, moves='birth,merge', sF=sf, doWriteStdOut=False)

    hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=100, nTask=1, K=K, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=sf, doWriteStdOut=False)

    observationModel = hmodel.obsModel
    numClusters = observationModel.K

    mus = [observationModel.get_mean_for_comp(k=i) for i in range(numClusters)]
    sigmas = [observationModel.get_covar_mat_for_comp(k=i) for i in range(numClusters)]

    target = bnpy.data.XData(X=(np.array(points)))
    LP = hmodel.calc_local_params(target)
    targetAssignments = np.argmax(LP['resp'], axis=1)

    LP = hmodel.calc_local_params(Data)
    fullAssignments = np.argmax(LP['resp'], axis=1)

    numPoints = []
    for i in range(numClusters):
		currX = np.array([Data.X[j] for j in range(len(Data.X)) if fullAssignments[j] == i])
		numPoints.append(currX.shape[0])

    return mus, sigmas, targetAssignments, numPoints, numClusters


def refineClustering(combo, assign, assignidx, samples, rdtol, baftol):
    assignment = {b : assign[assignidx[b]] for b in combo}
    clusters = set(assignment[b] for b in assignment)
    size = {c : float(sum(c == assignment[b] for b in combo)) for c in clusters}
    getbaf = (lambda c, p : float(sum(e[6] for b in combo for e in combo[b] if assignment[b] == c and e[0] == p)))
    baf = {c : {p : getbaf(c, p) / size[c] for p in samples} for c in clusters}
    getrdr = (lambda c, p : float(sum(e[1] for b in combo for e in combo[b] if assignment[b] == c and e[0] == p)))
    rdr = {c : {p : getrdr(c, p) / size[c] for p in samples} for c in clusters}

    mbaf = (lambda c : {p : baf[c][p] for p in samples})
    mrdr = (lambda c : {p : rdr[c][p] for p in samples})
    merge = {c : {'BAF' : mbaf(c), 'RDR' : mrdr(c), 'SIZE' : size[c], 'CLUS' : {c}} for c in clusters}

    def mergable(m):
        checkrdr = (lambda f, s : False not in set(abs(m[f]['RDR'][p] - m[s]['RDR'][p]) <= rdtol for p in samples))
        checkbaf = (lambda f, s : False not in set(abs(m[f]['BAF'][p] - m[s]['BAF'][p]) <= baftol for p in samples))
        check = (lambda f, s : checkrdr(f, s) and checkbaf(f, s))
        varrdr = (lambda f, s : sum(abs(m[f]['RDR'][p] - m[s]['RDR'][p]) for p in samples))
        varbaf = (lambda f, s : sum(abs(m[f]['BAF'][p] - m[s]['BAF'][p]) for p in samples))
        var = (lambda f, s : varrdr(f, s) + varbaf(f, s))

        seq = sorted(m, key=(lambda x : m[x]['SIZE']))
        for idx, f in enumerate(seq):
            opts = {s : var(f, s) for s in seq[idx+1:] for p in samples if s != f and check(f, s)}
            if len(opts) > 0:
                first = f
                break
        if len(opts) > 0:
            f = first
            s = sp.argmin(opts)
            return first, sp.argmin(opts)
        else:
            None

    m = mergable(merge)
    while m is not None:
        m1 = m[0]
        m2 = m[1]
        tot = float(merge[m1]['SIZE'] + merge[m2]['SIZE'])
        newbaf = {p : float(merge[m1]['BAF'][p] * merge[m1]['SIZE'] + merge[m2]['BAF'][p] * merge[m2]['SIZE']) / tot for p in samples}
        newrdr = {p : float(merge[m1]['RDR'][p] * merge[m1]['SIZE'] + merge[m2]['RDR'][p] * merge[m2]['SIZE']) / tot for p in samples}
        newclu = merge[m1]['CLUS'] | merge[m2]['CLUS']
        merge = {c : merge[c] for c in merge if c != m1 and c != m2}
        if len(merge) == 0:
            merge = {}
        merge[m1] = {'BAF' : newbaf, 'RDR' : newrdr, 'SIZE' : tot, 'CLUS' : newclu}
        m = mergable(merge)

    newassign = [-1 for i in range(len(assign))]
    for b in combo:
        get = [c for c in merge if assign[assignidx[b]] in merge[c]['CLUS']]
        assert len(get) == 1
        newassign[assignidx[b]] = get[0]
    assert -1 not in set(newassign)

    return newassign, len(merge)


def generateClouds(points, density, seed, sdeven=0.02, sdodd=0.02):
    res = []
    resappend = res.append
    for point in points:
        np.random.seed(seed=seed)
        for d in range(density):
            normal = np.random.normal
            newpoint = [normal(point[i], sdeven) if i%2==0 else normal(point[i], sdodd) for i in range(len(point))]
            resappend(newpoint)
    return res


def segmentBins(bb, clusters, samples):
    sbb = {bi : {record[0] : record[1:] for record in bb[bi]} for bi in bb}
    nbins = {cluster : {sample : len(clusters[cluster])  for sample in samples}  for cluster in clusters}
    rd = {cluster : {sample : float(sum(sbb[bi][sample][0] for bi in clusters[cluster])) / float(len(clusters[cluster])) for sample in samples} for cluster in clusters}
    nsnps = {cluster : {sample : sum(sbb[bi][sample][1] for bi in clusters[cluster]) for sample in samples} for cluster in clusters}
    cov = {cluster : {sample : float(sum(sbb[bi][sample][2] for bi in clusters[cluster])) / float(len(clusters[cluster])) for sample in samples} for cluster in clusters}
    return minSegmentBins(sbb, nbins, rd, nsnps, cov, clusters, samples)


def minSegmentBins(sbb, nbins, rd, nsnps, cov, clusters, samples):
    alpha = {cluster : {sample : sum(min(sbb[bi][sample][3], sbb[bi][sample][4]) for bi in clusters[cluster]) for sample in samples} for cluster in clusters}
    beta = {cluster : {sample : sum(max(sbb[bi][sample][3], sbb[bi][sample][4]) for bi in clusters[cluster]) for sample in samples} for cluster in clusters}
    mean = {cluster : {sample : float(alpha[cluster][sample]) / float(alpha[cluster][sample]+beta[cluster][sample]) for sample in samples} for cluster in clusters}
    return {cluster : {sample : (nbins[cluster][sample], rd[cluster][sample], nsnps[cluster][sample], cov[cluster][sample], alpha[cluster][sample], beta[cluster][sample], mean[cluster][sample]) for sample in samples} for cluster in clusters}


def scaleBAF(segments, samples, diploidbaf):
    diploid = -1
    main = -1
    for key in segments:
        if sum( (0.5 - segments[key][sample][-1]) <= diploidbaf for sample in samples) == len(samples):
            sample = list(samples)[0]
            if main < segments[key][sample][0]:
                main = segments[key][sample][0]
                diploid = key
    if diploid == -1:
        raise ValueError(sp.error("No potential neutral cluster has been found within the given threshold {}!".format(diploidbaf)))
    scalings = {sample : segments[diploid][sample][-1] for sample in samples}
    scale = (lambda value, scale, diploidbaf : min( (float(value) / float(scale)) * 0.5, 0.5) if (0.5-value) <= diploidbaf and scale > 0 else value)
    scaledBAF = {segment : {sample : scale(segments[segment][sample][-1], scalings[sample], diploidbaf) for sample in samples} for segment in segments}
    regulate = (lambda record, baf : record[:-3] + splitBAF(baf, record[-3] + record[-2]) + (baf,) if record[-1] != baf else record )
    return {segment : {sample : regulate(segments[segment][sample], scaledBAF[segment][sample]) for sample in samples} for segment in segments}


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


if __name__ == '__main__':
    main()
