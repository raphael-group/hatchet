#!/usr/bin/python3

import os, sys
import os.path
import argparse
import shlex
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value

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
    
    #if args["refpanel"] == "1000GP_Phase3":
       # wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz  
       # tar 1000GP_Phase3.tgz

if __name__ == '__main__':
    main()
