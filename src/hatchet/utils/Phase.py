#!/usr/bin/python3

import os, sys
import os.path
import argparse
import shlex
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
import requests
import tarfile

from . import ArgParsing as ap
from .Supporting import *
from . import Supporting as sp
from . import ProgressBar as pb



def main(args=None):
    log(msg="# log notes\n", level="STEP")
    args = ap.parse_phase_arguments(args)
    logArgs(args, 80)

    for i in args:
        print(i, args[i])

    if args["refpanel"] == "1000GP_Phase3":
        url = "https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz"
        out = args["outputphase"] + "/1000GP_Phase3.tgz"
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            f.write(r.content)
        # extract tar file, remove
        t = tarfile.open(out)
        t.extractall(args["outputphase"])
        t.close()
        os.remove(out)
        # tar 1000GP_Phase3.tgz

if __name__ == '__main__':
    main()
