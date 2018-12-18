#!/usr/bin/python2

import sys
import datetime
import os.path
import argparse
import shlex
import subprocess
from multiprocessing import Pool, Lock
import shutil


def parse_args():
    """
    Parse command line arguments
    Returns:
    """
    description = ""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-n","--normal", required=True, type=str, help="Indexed-sorted BAM file of mathed normal sample")
    parser.add_argument("-t","--tumors", required=True, type=str, nargs='+', help="A list of indexed-sorted BAM files of tumor samples")
    parser.add_argument("-s","--names", required=False, default=None, type=str, nargs='+', help="A list of sample names including normal as first and tumors in the same order as given to map the corresponding values in the total count and copy-number files")
    parser.add_argument("-p","--proportions", required=True, type=str, help="Tumor proportions in the mixture in the format...")
    parser.add_argument("-c","--copynumbers", required=True, type=str, help="An allele-specific copy-number list for all tumor clones with corresponding names")
    parser.add_argument("-u","--totalcounts", required=False, default=None, type=str, help="A list file expressing the total number of reads on each sample with corresponding names")
    parser.add_argument("-st","--samtools", required=False, default="", type=str, help="Path to \"samtools\" executable (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-T","--temp", required=False, default="./tmp", type=str, help="Temporary directory for the current execution (default: ./tmp)")
    parser.add_argument("-o","--output", required=False, default="./mixed_sample", type=str, help="Output name (default: ./mixed_sample.bam)")
    parser.add_argument("-q","--qualityreads", required=False, default=20, type=int, help="The minimum read quality to consider (default: 20)")
    parser.add_argument("-j","--jobs", required=False, default=1, type=int, help="Number of parallel jobs (default: 1)")
    parser.add_argument("-rs","--seed", required=False, default=12, type=int, help="Random seed (default: 12)")
    args = parser.parse_args()

    samtools = os.path.join(args.samtools, "samtools")
    if which(samtools) is None:
        raise ValueError(error("samtools has not been found or is not executable!"))

    if args.names != None and len(args.tumors) != len(args.names)+1:
        raise ValueError(error("The tumor names must be in the same number as tumor samples plus normal-sample name!"))

    if not os.path.isfile(args.normal):
        raise ValueError(error("The SNP file does not exist"))
    for t in args.tumors:
        if not os.path.isfile(t):
            raise ValueError(error("Tumor BAM file {} does not exist".format(t)))

    if not os.path.isfile(args.copynumbers):
        raise ValueError(error("Copy-number file {} does not exist".format(args.copynumbers)))
    if args.totalcounts != None and not os.path.isfile(args.totalcounts):
        raise ValueError(error("Totalcount file {} does not exist".format(args.totalcounts)))

    if os.path.exists(args.temp):
        shutil.rmtree(args.temp)
    os.makedirs(args.temp)

    if args.qualityreads <= 0:
        raise ValueError(error("The number of jobs must be a non-zero positive integer!"))
    if args.jobs <= 0:
        raise ValueError(error("The minium read quality must be a non-zero positive integer!"))
    if args.seed < 0:
        raise ValueError(error("The seed must be a non-zero positive integer!"))


    if not isinstance(args.tumors, list):
        tumors = [args.tumors]
    else:
        tumors = args.tumors
    if args.names == None:
        names = [os.path.splitext(os.path.basename(args.normal))[0]] + [os.path.splitext(os.path.basename(t))[0] for t in tumors]
    elif not isinstance(args.names, list):
        names = [args.names]
    else:
        names = args.names

    proportions = []
    for p in args.proportions.strip().split():
        sample =  [float(u) if u != '' else None for u in p.split(':')]
        if len(sample) != len(names):
            raise ValueError(error("{} is bad formatted!\nFor every tumor and normal sample the proportions of each clone must be specified separated by :, an empty entry is used when clone is not present!".format(p)))
        elif sum(e if e != None else 0 for e in sample) != 1.0:
            raise ValueError(error("The proportions of the following sample {} does not sum up to 1!".format(p)))
        else:
            proportions.append({name : u for name, u in zip(names, sample)})

    return {"normal" : args.normal,
            "tumors" : tumors,
            "names" : names,
            "proportions" : proportions,
            "totalcounts" : args.totalcounts,
            "copynumbers" : args.copynumbers,
            "samtools" : samtools,
            "temp" : args.temp,
            "output" : args.output,
            "jobs" : args.jobs,
            "q" : args.qualityreads,
            "seed" : args.seed}


def main():
    sys.stderr.write(log("# Parsing and checking arguments\n"))
    args = parse_args()
    sys.stderr.write(info("\n".join(["\033[96m{}:\t{}\033[0m".format(key, args[key]) for key in args]) + "\n"))

    sys.stderr.write(log("# Get total read counts for each sample\n"))
    counts = getTotalCounts(normal=args["normal"], tumors=args["tumors"], names=args["names"], samtools=args["samtools"], totalcounts=args["totalcounts"], j=args["jobs"], q=args["q"])
    for sample in counts: sys.stderr.write(info("In BAM of {} the total-read count is {}\n".format(sample, counts[sample])))
    
    sys.stderr.write(log("# Compute genome length for each clone\n"))
    lengths = computeLengths(copynumbers=args["copynumbers"], names=args["names"])
    for sample in lengths: sys.stderr.write(info("{} genome of length {}\n".format(sample, lengths[sample])))    

    for proportions in args["proportions"]:

        sys.stderr.write(log("# Generating bulk sample with proportions: {}\n".format(", ".join(["{}= {}".format(sample, proportions[sample]) for sample in proportions]))))

        sys.stderr.write(log("## Compute sampling proportions\n"))
        mixing = computeMixing(names=args["names"], counts=counts, lengths=lengths, proportions=proportions)
        for sample in mixing: sys.stderr.write(info("\033[96mFrom {} sampling {}\033[0m\n".format(sample, mixing[sample])))

        sys.stderr.write(log("## Sample from given BAM files\n"))
        samplings = sampling(normal=args["normal"], tumors=args["tumors"], names=args["names"], tmp=args["temp"], mixing=mixing, samtools=args["samtools"], q=args["q"], j=args["jobs"], seed=args["seed"])
        for sample, sampled in zip([args["normal"]]+args["tumors"], samplings): sys.stderr.write(info("Sampled {} from {}\n".format(sampled, sample)))

        output = os.path.splitext(os.path.basename(args["output"]))[0] + "_" + ("_".join(["{}{}".format(proportions[u], u) if u != None else 0 for u in proportions])).replace('.', '') + ".bam"
        sys.stderr.write(log("## Merging samples from BAM files in the resulting mixed BAM\n"))
        merge(output=output, normal=args["normal"], samtools=args["samtools"], samplings=samplings, j=args["jobs"], tmp=args["temp"])

        sys.stderr.write(log("## Removing temporary folder\n"))
        #        shutil.rmtree(tmp)
        
    return


def getTotalCounts(normal, tumors, names, samtools, totalcounts, j, q):
    counts = {}
    if totalcounts == None:
        pool = Pool(processes=min(j, len(names)))
        counts = pool.map(totalCount, [(samtools, name, sample, q) for (name, sample) in zip(names, [normal]+tumors)])
        counts = {name : count for (name, count) in counts}
    else:
        found = set()
        with open(totalcounts, 'r') as f:
            for line in f:
                if line != '':
                    line = line.strip().split()
                    found.add(line[0])
                    counts[line[0]] = int(line[1])
        if found != set(names):
            raise ValueError(error("The total read counts for the following sample names have not been found: {}".format(", ".join([name for name in names if not name in found]))))
    return counts


def totalCount(caps):
    samtools = caps[0]
    name = caps[1]
    sample = caps[2]
    quality = caps[3]
    cmd = "{} view {} -c -q {}".format(samtools, sample, quality)
    stdout, stderr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if stderr != "":
        raise ValueError(error("samtools warns \"{}\"on (sample={})".format(stderr, sample)))
    return (name, int(stdout.strip()))


def computeLengths(copynumbers, names):
    tumornames = names[1:]
    normalLength = 0
    lengths = {name : 0 for name in tumornames}
    idxs = {}
    with open(copynumbers, 'r') as f:
        for line in f:
            if line != '':
                parsed = line.strip().split()
                if line[0] == '#':
                    idxs = {clone : idx for idx, clone in enumerate(parsed) if idx > 2}
                    found = set(name for name in tumornames if name in idxs)
                    if(set(tumornames) != found):
                        raise ValueError(error("The following tumor names have not be founded in copy-number file: {}".format(", ".join([name for name in tumornames if not name in found]))))
                else:
                    normalLength += 2 * (int(parsed[2]) - int(parsed[1]))
                    lengths = {name : lengths[name] + (sum(map(int, parsed[idxs[name]].split('|'))) * (int(parsed[2]) - int(parsed[1]))) for name in tumornames}
    lengths[names[0]] = normalLength
    return lengths


def computeMixing(names, counts, lengths, proportions):
    target = min([counts[name] for name in counts])
    sys.stderr.write(info("The selected final total-read count is {}\n".format(target)))

    rescale = {sample : float(proportions[sample]) * float(lengths[sample]) for sample in names if proportions[sample] != None}
    norm = float(sum( float(rescale[sample]) for sample in rescale ))
    rescale = {sample : float(rescale[sample]) / float(norm) for sample in rescale}
    resize = {sample : float(target) / float(counts[sample]) for sample in rescale}

    return {sample : float(rescale[sample]) * float(resize[sample]) for sample in rescale}


def sampling(normal, tumors, names, tmp, mixing, samtools, q, j, seed):
    shutil.rmtree(tmp)
    os.makedirs(tmp)
    pj = min(j, len(mixing))
    pool = Pool(processes=pj)
    at = max(1, int( (j - pj) / len(mixing)))
    samples = {name : sample for name, sample in zip(names, [normal]+tumors) if name in mixing}
    samplings = pool.map(runsampling, [(samtools, samples[sample], tmp, mixing[sample], q, seed, at) for sample in samples])
    return samplings


def runsampling(box):
    samtools = box[0]
    sample = box[1]
    tmp = box[2]
    u = box[3]
    q = box[4]
    seed = box[5]
    at = box[6]
    name = os.path.splitext(os.path.basename(sample))[0] + "_sampled_{}.bam".format(str(u)[2:])
    name = os.path.join(tmp, name)
    cmd = "{} view {} -h -b -s {}.{} -q {} -o {} -@ {}".format(samtools, sample, seed, str(u)[2:], q, name, at)
    sys.stderr.write(info("Sampling {} from {}: {}\n".format(u, sample, cmd)))
    stdout, stderr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if stderr != "":
        raise ValueError("{}{}: samtools warns \"{}\"on (sample={}){}".format("\033[93m", self.name, stderr, sample, "\033[0m"))
    else:
        return name


def merge(output, normal, samtools, samplings, j, tmp):
    headerfile = os.path.join(tmp, "header.sam")
    cmd = "{} view -H {}".format(samtools, normal)
    sys.stderr.write(info("{}\n".format(cmd)))
    stdout, stderr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if stderr != "":
        sys.stderr.write(error("samtools warns \"{}\" on merging\n".format(stderr)))
    header = stdout + "@RG\tID:bulk_sample\tSM:None\tLB:None\tPL:Illumina\n"
    with open(headerfile, 'w') as f: f.write(header)

    cmd = "{} merge -@ {} -h {} -c -p {} {}".format(samtools, j, headerfile, output, " ".join(samplings))
    sys.stderr.write(info("{}\n".format(cmd)))
    stdout, stderr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if stderr != "":
        sys.stderr.write(error("samtools warns \"{}\" on merging\n".format(stderr)))

    sortedout = output.replace('.bam', '.sorted.bam')
    cmd = "{} sort -O bam -o {} -T {} -@ {} {}".format(samtools, sortedout, tmp, j, output)
    sys.stderr.write(info("{}\n".format(cmd)))
    stdout, stderr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if stderr != "":
        sys.stderr.write(error("samtools warns \"{}\" on sorting\n".format(stderr)))

    cmd = "{} index {}".format(samtools, sortedout)
    sys.stderr.write(info("{}\n".format(cmd)))
    stdout, stderr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if stderr != "":
        sys.stderr.write(error("samtools warns \"{}\" on indexing\n".format(stderr)))

    return


def error(msg):
    return "{}{}{}".format("\033[91m\033[1m", msg, "\033[0m")


def log(msg):
    timestamp = '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
    return "{}[{}]{}{}".format("\033[95m\033[1m", timestamp, msg, "\033[0m")


def info(msg):
    return "{}{}{}".format("\033[96m", msg, "\033[0m")


def which(program):
    import os

    def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file): return exe_file
    return None



if __name__ == '__main__':
        main()
