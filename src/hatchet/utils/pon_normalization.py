import re
import sys
import os
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy import optimize
from math import floor

def panel_normalize(bb, avg_normal, norm_constant):
    denoms = []
    panel_normalized = []
    normlens = []
    for _, r0 in bb.iterrows():
        bin_start = r0.START
        bin_end = r0.END
        bin_ch = r0['CHR']
        num = r0.TOTAL_READS

        # Find all panel bins overlapping this bin
        my_bins = avg_normal[(avg_normal.START <= bin_end) & (avg_normal.END >= bin_start) &
                            (avg_normal['CHR'] == bin_ch)]
        my_denom = 0
        my_bp = 0
        normlen = 0
        for _, r in my_bins.iterrows():
            if r.START < bin_start:
                if r.END < bin_end:
                    # suffix of panelbin overlaps bin
                    my_denom += (r.END - bin_start) * r.perbase
                    normlen += r.END - bin_start
                else:
                    # bin is within panelbin
                    my_denom += (bin_end - bin_start) + r.perbase
                    normlen += bin_end - bin_start

            else:
                if r.END > bin_end:
                    # prefix of panelbin overlaps bin
                    my_denom += (bin_end - r.START) * r.perbase
                    normlen += bin_end - r.START
                else:
                    # panelbin is within bin
                    my_denom += (r.END - r.START) * r.perbase
                    normlen += r.END - r.START
        normlens.append(normlen)

        denoms.append(my_denom)
        # if sample has more reads overall than PoN, downweight number of reads by norm_constant, where
        # norm_constant = (total reads in PoN) / (total reads in sample)
        
        try:
            panel_normalized.append( (num * norm_constant) / my_denom)
        except:
            import pdb; pdb.set_trace()

    return panel_normalized, normlens, denoms

def correct_baf(bb):
    for index, row in bb.iterrows():
        baf = row['BAF']
        coverage = floor(row['TOTAL_SNP_READS']/row['SNPS'])
        if coverage < 2:
            continue
        bb.at[index,'BAF'] = correct_baf_bin(baf,coverage)


def correct_baf_bin(baf, read_per_snp):
    # newton method functions f, f', and f''
    f = (lambda p: baf*(1-(1-p)**read_per_snp)-p)
    fprime = (lambda p: baf*read_per_snp*(1-p)**(read_per_snp-1)-1)
    fprime2 = (lambda p: -baf*read_per_snp*(read_per_snp-1)*(1-p)**(read_per_snp-2))
    initial_value=0.5
    return optimize.newton(f,initial_value,fprime=fprime,fprime2=fprime2)

def main(bb, outfile, pon_file):
    #bb_file = f"{outdir}/bb/bulk.bb"

    #bb = pd.read_csv(bb_file, sep='\t')
    all_normals = pd.read_csv(pon_file, sep="\t")
    pvt = all_normals.pivot(index = ['CHR', 'START'], columns = 'SAMPLE', values = 'RD').to_numpy()

    # normalize by the total read depth in each sample, then re-multiply to get comparable RDs
    total_per_sample = np.sum(pvt, axis = 0)
    avg_total = np.mean(total_per_sample) # axis = 0 sums along rows (e.g. takes sum of 1st col, then sum of 2nd col, ...)
    pvt = (avg_total * pvt) / total_per_sample

    # create a df with unique bins (PoN has 38 samples, 38 bins per locus)
    first_sample = all_normals.SAMPLE.unique()[0]
    avg_normal = all_normals[all_normals.SAMPLE == first_sample]
    avg_normal = avg_normal.sort_values(by = ['CHR', 'START'])
    avg_normal = avg_normal.drop(columns = 'SAMPLE')
    avg_normal.RD = all_normals.pivot(index = ['CHR', 'START'], columns = 'SAMPLE', values = 'RD').to_numpy().mean(axis = 1)

    norm_constant = avg_normal.RD.sum() / bb.TOTAL_READS.sum()
    avg_normal['CHR'] = avg_normal.CHR
    avg_normal['perbase'] = avg_normal.RD / (avg_normal.END - avg_normal.START)

    panel_normalized, normlens, denoms = panel_normalize(bb, avg_normal, norm_constant)

    #bb['normlen'] = normlens
    #bb['LENGTH'] = bb.END - bb.START
    bb = bb.drop(['RD'], axis=1)
    bb['RD'] = panel_normalized
    correct_baf(bb)
    bb['BAF'] = bb['BAF'].round(5)
    bb.to_csv(outfile, index=False, sep='\t',float_format="%.5f")
    #bb.to_csv(f"{outdir}/bb/bulk_renorm.bb", sep='\t', index=False)