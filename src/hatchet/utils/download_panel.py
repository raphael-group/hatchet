#!/usr/bin/python3

import os, sys
import os.path
import argparse
import subprocess as pr
import requests
import tarfile
import gzip
from . import ArgParsing as ap
from .Supporting import *
from . import Supporting as sp

def main(args=None):
    log(msg="# log notes\n", level="STEP")
    args = ap.parse_download_panel_arguments(args)
    logArgs(args, 80)

    if args["refvers"] not in ["hg19", "hg38"]:
        raise ValueError(sp.error("The reference genome version of your samples is not \"hg19\" or \"hg38\", please specify one of these two options!\n"))

    # download reference panel, prepare files for liftover
    if args["refpanel"] == "1000GP_Phase3":
        # download 1000GP ref panel
        if not os.path.isfile( os.path.join(args["refpaneldir"], "1000GP_Phase3", "1000GP_Phase3.sample") ):
            dwnld_1kgp(path = args["refpaneldir"])
        panel = os.path.join( args["refpaneldir"], "1000GP_Phase3" )
    else:
        raise ValueError(sp.error("Currently, only the 1000 genome panel aligned to GRCh37 without \"chr\" prefix is supported\n")) 
        #panel = args["refpanel"]       # for future: include custom panel

    # download necessary liftover files; 1000GP in hg19 coordinates
    # if users aligned reads to the hg38 build, we need to liftover coordinates to the reference panel (hg38 -> hg19)
    # since the 1000GP panel is in hg19 coordinates, we need to download (1) hg19  genome and (2) chain files
    # for liftover via picard
    hg19_path = dwnld_refpanel_genome(path=args["refpaneldir"])
    chains = dwnld_chains(path=args["refpaneldir"], refpanel_chr="False", sample_chr=args["chrnot"] )
    # if users aligned reads to the same reference genome as used in the reference panel, liftover isn't required, but
    # there could be different naming conventions of chromosomes, with or without the 'chr' prefix. The 1000GP reference 
    # panel does not use 'chr' prefix, so input into shapeit also should not have this
    rename_files = mk_rename_file(path = args["refpaneldir"]) 

def dwnld_chains(path, refpanel_chr, sample_chr):
    # order of urls important! [0] chain for sample -> ref panel, [1] chain for ref panel -> sample
    urls = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
    paths = [os.path.join(path, os.path.basename(i)) for i in urls] # paths for url downloads
    for i, url in enumerate(urls):
        r = requests.get(url, allow_redirects=True)
        open(paths[i], 'wb').write(r.content)
    # make all necessary chain files to convert from hg38 (w/ or w/out chr notation) to hg19 (no chr notation),
    # and also to lift back over from hg19 (no chr notation) to hg38 (w/ or w/out chr notation).
    c1_1 = mod_chain(paths[0], refpanel_chr, sample_chr="True", refpanel_index=7, sample_index=2) # modify chr notation of hg38ToHg19, ref panel chr in 7th field, sample chr in 2nd field
    c1_2 = mod_chain(paths[0], refpanel_chr, sample_chr="False", refpanel_index=7, sample_index=2) 
    c2_1 = mod_chain(paths[1], refpanel_chr, sample_chr="True", refpanel_index=2, sample_index=7) # modify chr notation of hg19ToHg38, ref panel chr in 2th field, sample chr in 7nd field
    c2_2 = mod_chain(paths[1], refpanel_chr, sample_chr="False", refpanel_index=2, sample_index=7) 
    [ os.remove(p) for p in paths ] # remove original chain files
    return 

def mod_chain(infile, refpanel_chr, sample_chr, refpanel_index, sample_index):
    if sample_chr == "True":
        name = infile.strip(".gz").replace("over", "chr")
    elif sample_chr == "False":
        name = infile.strip(".gz").replace("over", "no_chr")

    with open(name, 'w') as new:
        with gzip.open(infile, 'rt') as f:
            for l in f:
                if l.startswith("chain"):
                    l = l.split()
                    if refpanel_chr == "False": l[refpanel_index] = l[refpanel_index].replace("chr","") 
                    if sample_chr == "False": l[sample_index] = l[sample_index].replace("chr","") 
                    new.write(" ".join(l) + "\n")
                else:
                    new.write(l)
    return name

def dwnld_refpanel_genome(path):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    out = os.path.join(path, "hg19.fa.gz")
    newref = os.path.join(path, "hg19_no_chr.fa")
    if not os.path.isfile(newref):
        # download hg19
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out, 'wb') as f:                                                                                                                    
                for chunk in r.iter_content(chunk_size=8192):                                                                                             
                    f.write(chunk)      
        # change chr notation
        with open(newref, 'w') as new:
            with gzip.open(out, 'rt') as f:
                [ new.write(l.replace("chr","")) if l.startswith(">") else new.write(l) for l in f ]
        os.remove(out)
        # make dict file
        dict_file = newref.replace(".fa",".dict")
        cmd = f"samtools dict {newref} > {dict_file}"
        errname = os.path.join(path, f"samtools.log")
        with open(errname, 'w') as err:
            run = pr.run(cmd, stdout=err, stderr=err, shell=True, universal_newlines=True)
        if run.returncode != 0:
            raise ValueError(sp.error(f"Samtools dict creation failed, please check errors in {errname}!"))
        else:
            os.remove(errname)
    return newref

def mk_rename_file(path):
    # makes rename_chrs1.txt for removing "chr", rename_chrs2.txt for adding "chr"
    names = [ os.path.join(path, f"rename_chrs{i}.txt") for i in range(1,3) ]
    for i, n in enumerate(names):
        with open(n, 'w') as f:
            for j in range(1,23):
                f.write(f"chr{j} {j}\n") if i == 0 else f.write(f"{j} chr{j}\n") 
    return names

def dwnld_1kgp(path):
    url = "https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz"
    out = os.path.join(path + "1000GP_Phase3.tgz")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(out, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        #f.write(r.content)
    # extract tar file, remove
    t = tarfile.open(out)
    t.extractall(path)
    t.close()
    os.remove(out)

if __name__ == '__main__':
    main()
