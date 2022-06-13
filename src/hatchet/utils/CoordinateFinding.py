#!/usr/bin/python


from math import floor
from math import ceil


lhg19 = [
    -1,
    249904550,
    243199373,
    198022430,
    191535534,
    180915260,
    171115067,
    159321559,
    146440111,
    141696573,
    135534747,
    135046619,
    133851895,
    115169878,
    107349540,
    102531392,
    90354753,
    81529607,
    78081510,
    59380841,
    63025520,
    48157577,
    51304566,
]


def extractChr(ref):
    return int(''.join([i for i in ref if i.isdigit()]))


def findStart(bamfile, seq, least=0):
    L = 0
    R = lhg19[extractChr(seq)]
    first = True

    while R - L > 1:
        if first:
            M = int(ceil(R * 0.01))
        else:
            M = int(ceil((L + R) / 2))
        count = bamfile.count(seq, L, M)
        # print "[{}, {}] = {}".format(L, M, count)
        if count > least:
            R = M
        else:
            L = M
        first = False

    if bamfile.count(seq, L, L) > least:
        return L
    else:
        return R


def findEnd(bamfile, seq, least=0):
    L = 0
    R = lhg19[extractChr(seq)]
    first = True

    while R - L > 1:
        if first:
            M = int(floor(R * 0.99))
        else:
            M = int(floor((L + R) / 2))
        count = bamfile.count(seq, M, R)
        # print "[{}, {}] = {}".format(M, R, count)
        if count > least:
            L = M
        else:
            R = M
        first = False

    if bamfile.count(seq, R, R) > least:
        return R
    else:
        return L


def binChr(bamfile, sample, seq, size, start=0, end=0, least=-1):
    res = []
    index = start
    if end == 0:
        end = lhg19[extractChr(seq)]
    while index < end:
        border = min(end, index + size)
        count = bamfile.count(seq, index, border)
        if count > least:
            res.append(
                {
                    'sample': sample,
                    'chromosome': extractChr(seq),
                    'start': index,
                    'end': (index + size),
                    'count': count,
                }
            )
        index += size
    return res
