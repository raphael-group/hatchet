#!/usr/bin/python


import sys, os, argparse

def parse_arguments():
    description = "Rescale a sample and encode in the mixcnp format for a single sample"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-s","--samples", required=True, nargs='+', type=str, help="The filenames of rescaled samples")
    args = parser.parse_args()
    filenames = args.samples
    return filenames


def main():
    samplesfiles = parse_arguments()
    samples = []
    starts = []
    ends = []

    for filename in samplesfiles:
        with open(filename) as f:
            lines = f.readlines()
            starts.append([map(int, i.split()) for i in lines[0].split('|')])
            ends.append([map(int, i.split()) for i in lines[1].split('|')])
            samples.append([map(float, i.split()) for i in lines[2].split('|')])

    num_samples = len(samples)
    if (len(starts) != num_samples or len(ends) != num_samples):
        raise Exception("Wrong number of samples between samples, starts, and ends!")

    num_chro = len(samples[0])

    for sample in range(num_samples):
        for chro in range(num_chro):
            previous = -1
            if(len(starts[sample][chro]) != len(ends[sample][chro]) or len(starts[sample][chro]) != len(samples[sample][chro])):
                raise Exception("Input format is wrong!")
            for seg in range(len(samples[sample][chro])):
                if(not(starts[sample][chro][seg] <= ends[sample][chro][seg])):
                    raise Exception("Found a segment's start that is following his end!")
                if(not(starts[sample][chro][seg] >= previous)):
                    raise Exception("The segments are not consecutive and disjoint")
                previous = ends[sample][chro][seg]

    breakpoints = [[] for i in range(num_chro)]

    for sample in range(num_samples):
        if(len(starts[sample]) != num_chro or len(ends[sample]) != num_chro):
            raise Exception("Found sample starts with wrong number of chromosomes")
        for chro in range(num_chro):
            for i in range(len(starts[sample][chro])):
                if starts[sample][chro][i] not in breakpoints[chro]:
                    breakpoints[chro].append(starts[sample][chro][i])
                if ends[sample][chro][i] not in breakpoints[chro]:
                    breakpoints[chro].append(ends[sample][chro][i])

    for bk in breakpoints:
        bk.sort()

    resulting_starts = []
    resulting_ends = []

    for chro in range(num_chro):
        temp_starts = []
        temp_ends = []
        for bk in range(1, len(breakpoints[chro])):
            start = breakpoints[chro][bk - 1]
            end = breakpoints[chro][bk]
            flag_global = True
            for sample in range(num_samples):
                flag_sample = False
                for seg in range(len(samples[sample][chro])):
                    if(start >= starts[sample][chro][seg] and end <= ends[sample][chro][seg]):
                        flag_sample = True
                if(not flag_sample):
                    flag_global = False
            if(flag_global):
                temp_starts.append(start)
                temp_ends.append(end)
        resulting_starts.append(temp_starts)
        resulting_ends.append(temp_ends)

    resulting_segments = [len(chro) for chro in resulting_starts]
    resulting_samples = []
    for sample in range(num_samples):
        temp_sample = []
        for chro in range(num_chro):
            temp_chro = []
            for seg in range(resulting_segments[chro]):
                temp_seg = -1
                for old in range(len(samples[sample][chro])):
                    if(resulting_starts[chro][seg] >= starts[sample][chro][old] and resulting_ends[chro][seg] <= ends[sample][chro][old]):
                        temp_seg = samples[sample][chro][old]
                temp_chro.append(temp_seg)
            temp_sample.append(temp_chro)
        resulting_samples.append(temp_sample)

    sys.stderr.write(" | ".join([' '.join(map(str, i)) for i in resulting_starts]) + "\n")
    sys.stderr.write(" | ".join([' '.join(map(str, i)) for i in resulting_ends]) + "\n")

    print("#PARAMS")
    print(str(num_chro) + " #number of chromosomes")
    print(str(num_samples) + " #number of samples")
    print(" ".join(map(str, resulting_segments)) + " #number of segments for each chromosome")
    print("#SAMPLES")
    counter = 0
    for sample in resulting_samples:
        print(str(counter) + " : " + " | ".join([' '.join(map(str, i)) for i in sample]))
        counter += 1


if __name__ == '__main__':
	main()
