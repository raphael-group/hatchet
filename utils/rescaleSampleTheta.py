#!/usr/bin/python


import sys, os, argparse
import pandas as pd


def parse_arguments():
    description = "Rescale a sample and encode in the mixcnp format for a single sample"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-s","--sample", required=True, type=str, help="The filename of a sample in Theta's format")
    parser.add_argument("-S","--rescalingSegments", nargs='+', required=True, type=str, help="The id of diploid segments used for rescaling")
    args = parser.parse_args()
    filename = args.sample
    rescalingSegments = args.rescalingSegments
    return filename, rescalingSegments


def main():
    filename, rescaling = parse_arguments()

    data = []
    factor = 0
    norm = 0
    counter = 0
    fractions = ""
    starts = ""
    ends = ""
    with open(filename) as f:
        for l in f:
            line = l.strip().split()
            data.append({"id" : line[0], "chr" : int(line[1]), "start" : int(line[2]), "end" : int(line[3]), "tumorCount" : int(line[4]), "normalCount" : int(line[5])})

    df = pd.DataFrame(data)
    df.sort_values(by=['chr', 'start'], ascending=[True, True])
    df["ratio"] = df["tumorCount"] / df["normalCount"]

    for index, row in df.iterrows():
        if(not(row["start"] <= row["end"])):
            raise Exception("A segment's start is coming after his end")
        if row["id"] in rescaling:
            factor += row["ratio"]
            counter += 1
    factor = float(factor / counter)

    if(counter != len(rescaling)):
        raise Exception("Not all the diploid rescaling segments have been found!")

    for chro in range(1,23):
        slicer = df[df["chr"]==chro]
        previous = -1
        for index, row in slicer.iterrows():
            f = float(2 * row["tumorCount"]) / float(factor * row["normalCount"])
            starts += str(row["start"]) + " "
            ends += str(row["end"]) + " "
            fractions += str(f) + " "
            if(not(row["start"] >= previous)):
                raise Exception("The segments are not consecutive and disjoint in the same sample!")
            previous = row["end"]
        if(chro < 22):
            starts += "| "
            ends += "| "
            fractions += "| "

    print(starts)
    print(ends)
    print(fractions)


if __name__ == '__main__':
	main()
