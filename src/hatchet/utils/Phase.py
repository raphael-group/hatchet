#!/usr/bin/python3

import os, sys
import os.path
import argparse
import shlex
import subprocess as pr
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
from scipy.stats import beta
# from statsmodels.stats.proportion import *

from . import ArgParsing as ap
from .Supporting import *
from . import Supporting as sp
from . import ProgressBar as pb



def main(args=None):
    log(msg="# log notes\n", level="STEP")

    print("hello world")

if __name__ == '__main__':
    main()
