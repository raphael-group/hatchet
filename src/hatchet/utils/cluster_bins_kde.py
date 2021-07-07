#!/usr/bin/python3

import os, shutil
from os.path import expanduser
import sys
import math
import copy
import numpy as np
from scipy.stats import multivariate_normal, poisson, beta, binom, gaussian_kde
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict, Counter

from .ArgParsing import parse_cluster_kde_args
from . import Supporting as sp


def main(args=None):
    sp.log(msg="# Parsing and checking input arguments\n", level="STEP")
    # TODO: implement in ArgParsing
    args = parse_cluster_kde_args(args)
    
    sp.log(msg="# Reading the combined BB file\n", level="STEP")
    
    bb = pd.read_table(args["bbfile"])
    if 'CHR' in bb:
        # required for compatibility with
        raise IndexError(sp.error("BB file chromosome column name should be '#CHR', not 'CHR'."))
        #chr_col = 'CHR'
    else:
        #chr_col = '#CHR'
        pass
    
    chr_col = '#CHR'
    if '#CHR' not in bb:
        raise IndexError(sp.error("BB file is missing '#CHR' column!"))
    
    bb = bb.sort_values(['SAMPLE', chr_col, 'START']) 
    bb = bb.reset_index(drop = True)
    if len(bb.SAMPLE.unique()) > 1:
        raise ValueError(sp.error("ERROR: This dataframe has multiple samples, but cluBB_KDE supports only single-sample data."))

    required_columns = ['RD', 'ALPHA', 'BETA', 'TOTAL_READS']
    missing_columns = [a for a in required_columns if a not in bb]
    if len(missing_columns) > 0:
        raise ValueError(sp.error("Missing columns from BB that are required for modeling: {}".format(missing_columns)))


    sp.log(msg="# Removing outlier bins\n", level="STEP")
    cap = 2 * np.percentile(bb.RD, 99.9)
    outliers = np.where(bb.RD > cap)[0]
    bb_kde = bb.drop(index = outliers)
    #bb = bb.reset_index(drop = True)
    sp.log(msg=f"Removed {len(outliers)} outlier bins for KDE centroids\n", level = 'INFO')

    sp.log(msg="# Identifying centroids using KDE\n", level="STEP")
    if "UNCORRECTED" in bb_kde:
        arr = np.array([np.abs(0.5 - bb_kde.UNCORRECTED), bb_kde.RD]).T    
    else:
        arr = np.array([np.abs(0.5 - bb_kde.BAF), bb_kde.RD]).T
    _, centers, _, _, _ = kde_centers_gridfit(arr, min_center_density = args['centroiddensity'], 
                                              bandwidth = args['bandwidth'], grid_dim = args['mesh'], 
                                              yvar = args['variance'], min_grid_density = args['griddensity'],
                                              max_copies = args['maxcopies'], fname = args['outfigure'], verbose = True)

    


    sp.log(msg="# Assigning bins to centroids\n", level="STEP")
    if args['snpsfile'] is not None:
        snps = read_SNPs_iter(args['snpsfile'], bb)
    
    # Re-assign until all clusters are at least size [-s, --minsize]
    min_cluster_size = args['minsize']
    while True:
        if args['snpsfile'] is not None:
            labels, affinities, _, _, bad_bins = assign_bins_binom(bb, snps, centers, chr_col)
            if len(bad_bins) > 0:
                raise ValueError(sp.error(f"Found {len(bad_bins)} bins with no SNPs: {bad_bins}"))
        else:
            labels, affinities, _, _ = assign_bins_beta(bb, centers)
    
        cntr = Counter(labels)
        if cntr.most_common()[-1][1] >= min_cluster_size:
            break
        else:
            keep_ids = [i for (i, v) in Counter(labels).items() if v >= min_cluster_size]
            centers = centers[keep_ids]

    # Reindex labels to be numbered from 1 to n_clusters
    labels = reindex(labels)

    # Form bbc output file
    columns = [chr_col, 'START', 'END', 'SAMPLE', 'RD', '#SNPS', 'COV', 'ALPHA', 'BETA', 'BAF']
    bbc = bb.copy()[columns]
    bbc = bbc.rename(columns=lambda x: '#CHR' if x == 'CHR' else x)
    bbc['CLUSTER'] = labels
    bbc.to_csv(args['outbins'], index = False, sep = '\t')

    # Form seg output file
    cluster_dfs = [df for df in bbc.groupby("CLUSTER")]
    rows = []
    d = args['diploidbaf']
    for i, df in cluster_dfs:
        assert len(df.SAMPLE.unique()) == 1
        row = [i, df.SAMPLE.unique()[0], len(df), df.RD.mean(), df['#SNPS'].sum(), df.COV.mean(), 
               df.ALPHA.sum(), df.BETA.sum(), df.BAF.mean()]
        if row[-1] > 0.5 - d:
            row[-1] = 0.5
        rows.append(row)        
    seg = pd.DataFrame(rows, columns = ['#ID', 'SAMPLE', '#BINS', 'RD', '#SNPS', 'COV', 'ALPHA', 'BETA', 'BAF'])
    seg.to_csv(args['outsegments'], index = False, sep = '\t')

def read_SNPs_iter(snpsfile, bb, suppress_warnings = True):  
    """
    Assumes that bins in both the snps dataframe and bins in bb are sorted in the same order
    """
    
    snps = pd.read_csv(snpsfile, delimiter = '\t', 
                       header = None, names = ['CHR', 'POS', 'SAMPLE', 'REF', 'ALT'],
                       dtype = {'CHR':np.object, 'POS':np.uint64, 'SAMPLE':np.object, 'REF':np.uint32, 'ALT':np.uint32}).sort_values(['CHR', 'POS'])     
    snps = snps.reset_index(drop = True)
    
    # sorted list of bins for each chromosome
    bins = {str(ch):sorted(set([(row.START, row.END) for _, row in bb[bb['#CHR'] == ch].iterrows()])) for ch in bb['#CHR'].unique()}
    
    # mapping from chromosome to current index in the list of bins
    chr2pos = {ch:0 for ch in snps.CHR.unique()}
    
    data = defaultdict(list)
    
    sp.log(msg= f"Reading SNPs file: {snpsfile}\n", level = "INFO")
    for i, row in snps.iterrows():
        # TODO: progressbar that moves in increments of 100k rows
                
        sample = row.SAMPLE
        ch = row.CHR
        pos = row.POS
        ref = row.REF
        alt = row.ALT
        
        if ch not in bins:
            if not (ch.endswith('X') or ch.endswith('Y')):
                sp.log(msg=f"WARNING: missing chromosome {ch}.\n", level = "WARN")
            continue
    
        # find the bin that contains this SNP
        j = chr2pos[ch]
        while j < len(bins[ch]) and pos > bins[ch][j][1]:
            j += 1
        
        if j >= len(bins[ch]):
            if not suppress_warnings:
                sp.log(msg=f"WARNING: CHR {ch}, POS {pos} is beyond bins (perhaps wrong centromeres file).\n", level = "WARN")
            
            continue
            
        start = bins[ch][j][0]
        end = bins[ch][j][1]
        if pos < start or pos > end:
            if not suppress_warnings:
                sp.log(msg=f"WARNING: CHR {ch}, POS {pos} is not in any bin (perhaps wrong centromeres file)\n", level = "WARN")
            continue
            #print(pos, start, end, ch, j)
            #return bins[ch]
        
        chr2pos[ch] = j


        data[ch, start, end].append((sample, pos, ref, alt))

    
    sp.log(msg="Done reading SNPs file, constructing dataframes...\n", level = "INFO")
    
    result = {}
    for b, slist in data.items():
        if len(slist) == 0:
            print(b)
            continue
        sample, pos, ref, alt = zip(*slist)
        df = pd.DataFrame()
        df['SAMPLE'] = sample
        df['POS'] = pos
        df['ALT'] = alt
        df['REF'] = ref
        df['TOTAL'] = np.array(alt) + np.array(ref)
        
        result[b] = df

    sp.log("Finished constructing dataframe\n", level = "INFO")

    return result

def assign_bins_beta(bins, centers):
    """
    Assign bins according to a probabilistic model, using a beta distribution for BAFs.
    """
    # Bins is a dataframe with columns ID, READS, NORMAL_READS,
    # SNPs is a dictionary that maps bin ID to a dataframe with columns ALPHA, BETA w/ alternate and reference read counts
    
    # Assign each bin to the cluster with maximum probability
    affinities = []
    term1s = []
    term2s = []
    
    for _, row in bins.iterrows():        
        t = row.CORRECTED_READS
        y = row.NORMAL_READS

        a = row.ALPHA
        b = row.BETA
        
        term1 = poisson.logpmf(t, y * centers[:,1])
        term2 = beta.logpdf(0.5 - centers[:, 0], a, b)
        
        affinities.append(term1 + term2)
        term1s.append(term1)
        term2s.append(term2)

    
    affinities = np.array(affinities)
    labels = np.argmax(affinities, axis = 1)
    return labels, affinities, term1s, term2s

def assign_bins_binom(bins, snps, centers, chr_col):
    """
    Assign bins according to a probabilistic model, using a binomial distribution for SNP read counts.
    """
    # Bins is a dataframe with columns ID, READS, NORMAL_READS, and LENGTH
    # SNPs is a dictionary that maps bin ID to a dataframe with columns ALT, REF w/ alternate and reference read counts
    
    
    # Assign each bin to the cluster with maximum probability
    affinities = []
    bad_bins = []
    term1s = []
    term2s = []
    
    i = 0
    for _, row in bins.iterrows():
        t = row.CORRECTED_READS
        y = row.NORMAL_READS
        bin_id = (row[chr_col], row.START, row.END)
        
        if bin_id not in snps:
            bad_bins.append(bin_id)
            affinities.append([0] * len(centers))
            continue
        
        alts = np.minimum(snps[bin_id].ALT, snps[bin_id].REF)
        totals = snps[bin_id].REF + snps[bin_id].ALT
        term1 = poisson.logpmf(t, y * centers[:,1])
        


        term2 = binom.logpmf(np.tile(alts, len(centers)), np.tile(totals, len(centers)), 
                             np.repeat(0.5 - centers[:,0], len(alts))).reshape(len(centers), len(alts)).sum(axis = 1)
        affinities.append(term1 + term2)
        term1s.append(term1)
        term2s.append(term2)

    
    affinities = np.array(affinities)
    labels = np.argmax(affinities, axis = 1)
    return labels, affinities, term1s, term2s, bad_bins

def kde_centers_gridfit(arr, min_center_density, bandwidth, grid_dim, yvar, 
                        min_grid_density, max_copies, fname, verbose):
    x = arr[:, 0]
    y = arr[:, 1]
 
    xmin, ymin = np.min(arr, axis = 0)
    ymax = np.max(y)
    #_, ymax = np.max(arr, axis = 0)
    xmax = 0.5
    
    # Compute KDE
    if verbose:
        sp.log(msg="Computing KDE\n", level = "INFO")
    xy = np.vstack([x,y])
    if bandwidth:
        kde = gaussian_kde(xy, bw_method = bandwidth)
    else:
        kde = gaussian_kde(xy)
    z = kde(xy)

    # Construct grid
    xvals = np.linspace(xmin, xmax, num = grid_dim)
    yvals = np.linspace(ymin, ymax, num = grid_dim)
    points = np.array(np.meshgrid(xvals, yvals)).T.reshape(-1, 2)

    # Evaluate KDE at each grid vertex
    if verbose:
        sp.log(msg="Evaluating grid\n", level = "INFO")
    probs = kde(points.T)
    probs_grid = probs.reshape(grid_dim, grid_dim)
    
    
    ## Identify local maxima ##
    lm = peak_local_max(probs_grid, footprint = np.ones((3,3)))
    lm_orig = lm.copy()
    
    if min_center_density == 'mean':
        min_center_density = np.mean(probs_grid[lm])
    
    ############### Grid analysis ########################
    _, purity, scaling, _ = find_grid(arr, variance = yvar)
    means = compute_means(purity, scaling, max_copies = max_copies)
    
    xvar = yvar / ( (np.max(y) - np.min(y)) /  (np.max(x) - np.min(x)) )

    vertices = vertex_gaussians(means, xvar, yvar)

    # Get coordinates of local optima
    candidates = np.array([xvals[lm[:, 0]], yvals[lm[:, 1]]]).T

    # For each vertex, find the index of the best optimum
    resp = np.array([vertices[i].pdf(candidates) for i in range(len(vertices))]).T
    best_idx = np.argmax(resp, axis = 0)
    best_probs = resp[best_idx, np.arange(resp.shape[1])]
        
    # Keep those best matches that have sufficient density (i.e., fit the grid sufficiently well)
    grid_idx = best_idx[np.where(best_probs > min_grid_density)[0]]
 
    # Map selected optima back to grid indices (to match lm)
    lm_grid = np.array([lm[grid_idx, 0], lm[grid_idx, 1]]).T
    
    # Remove local maxima below a certain density
    lm = lm[np.where(probs_grid[lm[:, 0], lm[:, 1]] > min_center_density)]
    
    # Take the union of the two sets as the final solution
    result = np.concatenate([lm, lm_grid])
    result = np.unique(result, axis = 0)
    
    if verbose:
        sp.log(msg=f"Found optima. Original: {len(lm_orig)}, grid: {len(lm_grid)}, threshold: {len(lm)}, union: {len(result)}\n", level = "INFO")
    
    if fname is not None:
        ### Mark local maxima on grids for visualization
        # Original local maxima (no thresholding)
        probs_show_orig = probs_grid.copy()
        probs_show_orig[lm_orig[:, 0], lm_orig[:, 1]] = np.max(probs_grid)

        # Local maxima that are also above density threshold
        probs_show_threshold = probs_grid.copy()
        probs_show_threshold[lm[:, 0], lm[:, 1]] = np.max(probs_grid)

        # Local maxima that are also close to grid vertices
        probs_show_grid = probs_grid.copy()
        probs_show_grid[lm_grid[:, 0], lm_grid[:, 1]] = np.max(probs_grid)

        # Local maxima present in the union
        probs_show_union = probs_grid.copy()
        probs_show_union[lm[:, 0], lm[:, 1]] = np.max(probs_grid)
        probs_show_union[lm_grid[:, 0], lm_grid[:, 1]] = np.max(probs_grid)

        plt.figure(dpi = 300)
        plt.subplot(221)
        plt.imshow(probs_show_orig.T[::-1])
        plt.title("Original local maxima")
        plt.subplot(222)
        plt.title("After density threshold")
        plt.imshow(probs_show_threshold.T[::-1])
        plt.subplot(223)
        plt.title("After grid threshold")
        plt.imshow(probs_show_grid.T[::-1])
        plt.subplot(224)
        plt.title("Union")
        plt.imshow(probs_show_union.T[::-1])
        plt.tight_layout()
        plt.savefig(fname)


    centers_x = xvals[result[:, 0]]
    centers_y = yvals[result[:, 1]]
    
    # Return the original coordinates as centers
    if verbose:
        sp.log(msg="Done finding centroids\n", level = "INFO")

    return z, np.array([centers_x, centers_y]).T, lm_orig, lm, probs_grid

def find_grid(arr, variance):     
    if np.median(arr[:,0]) > 0.25:
        sp.log("Mirroring BAF in find_grid\n", level = "INFO")
        # input is probably straight BAF instead of allelic imbalance
        # convert to allelic imbalance
        arr[:, 0] = np.abs(0.5 - arr[:, 0])
    rdr = arr[:, 1]
    imbalance = arr[:, 0]
    
    # Possible purity values 
    purities = np.linspace(0.05, 0.975, 38)

    x_variance = variance / (np.std(rdr) /  np.std(imbalance))
    
    # Copy number states for grid vertices
    cns = [(a, b) for a in range(1, 5) for b in range(a + 1)]
        
    # lambda functions from Simone's code
    # p is purity, float in (0, 1]
    # cn is a tuple of allele-specific copy numbers (A,B)
    # s is the scaling factor 
    calcFraction = (lambda p, cn : float(2.0 * (1.0 - p) + sum(cn) * p))
    calcRDR = (lambda p, cn, s : calcFraction(p, cn) / float(s) )
    calcBAF = (lambda p, cn : float(1.0 * (1.0 - p) + min(cn) * p) / calcFraction(p, cn) )   
    
    # not sure what d is
    #alcScalingFactor = (lambda p, d : float(2.0 + 2.0 * p) / float(d))
    
    best_score = -math.inf
    best_params = None
    idx = np.array(range(len(rdr)))
    
    for purity in purities:
        # Compute possible scalings using various RDR values 
        #scalings = [calcScalingFactor(purity, r) for r in rdr_candidates]
        scalings = np.linspace(0.1, 2.5, 121)
        
        for scaling in scalings:
            # Compute the locations of grid vertices
            means = np.array([(0.5 - calcBAF(purity, cn), calcRDR(purity, cn, scaling)) for cn in cns])
            
            # Put a diagonal Gaussian on each grid vertex
            vs = vertex_gaussians(means, x_variance, variance) 
            # Compute responsibilities            
            resp = np.array([v.logpdf(arr) for v in vs]).T
            
            # Compute log-likelihood of data
            labels = np.argmax(resp, axis = 1)
            
            #print(resp.shape, labels.shape)

            score = np.sum(resp[idx, labels])
            
            if score > best_score:
                if not best_params is None:
                    #print(best_score, best_params[1:3])
                    pass
                best_score = score
                best_params = means, purity, scaling, labels
            
    return best_params

def compute_means(purity, scaling, max_copies):
    cns = [(a, b) for a in range(1, max_copies) for b in range(a + 1)]
        
    calcFraction = (lambda p, cn : float(2.0 * (1.0 - p) + sum(cn) * p))
    calcRDR = (lambda p, cn, s : calcFraction(p, cn) / float(s) )
    calcBAF = (lambda p, cn : float(1.0 * (1.0 - p) + min(cn) * p) / calcFraction(p, cn) )   
    calcScalingFactor = (lambda p, d : float(2.0 + 2.0 * p) / float(d))

    return np.array([(0.5 - calcBAF(purity, cn), calcRDR(purity, cn, scaling)) for cn in cns])

def vertex_gaussians(means, x_variance, y_variance):
    return [multivariate_normal(mean = means[i], cov = [[x_variance, 0], [0, y_variance]]) for i in range((len(means)))]

def reindex(labels):
    """
    Given a list of labels, reindex them as integers from 1 to n_labels
    Labels are in nonincreasing order of prevalence
    """
    old2new = {}
    j = 1
    for i, _ in Counter(labels).most_common():
        old2new[i] = j
        j += 1
    old2newf = lambda x: old2new[x]

    return [old2newf(a) for a in labels]

if __name__ == '__main__':
    main()
