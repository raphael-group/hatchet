import sys
import datetime
import requests


def naturalOrder(text):
    return [int(s) if s.isdigit() else ord(s) for s in text]


def numericOrder(text):
    return int(digits(text))


def digits(string):
    return "".join(x for x in string if x.isdigit())


def argmax(d):
    return max(d, key=(lambda x : d[x]))


def argmin(d):
    return min(d, key=(lambda x : d[x]))


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
            if is_exe(exe_file):
                return exe_file

    return None


def urlexists(path):
    r = requests.head(path)
    return r.status_code == requests.codes.ok


def log(msg, level=None, lock=None):
    timestamp = '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
    if level == "STEP":
        if lock is None:
            sys.stderr.write("{}{}[{}]{}{}".format(bcolors.BOLD, bcolors.HEADER, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}{}[{}]{}{}".format(bcolors.BOLD, bcolors.HEADER, timestamp, msg, bcolors.ENDC))
    elif level == "INFO":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.OKGREEN, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.OKGREEN, timestamp, msg, bcolors.ENDC))
    elif level == "WARN":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.WARNING, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.WARNING, timestamp, msg, bcolors.ENDC))
    elif level == "PROGRESS":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.BBLUE, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.BBLUE, timestamp, msg, bcolors.ENDC))
    elif level == "ERROR":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.FAIL, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.FAIL, timestamp, msg, bcolors.ENDC))
    else:
        if lock is None:
            sys.stderr.write("{}".format(msg))
        else:
            with lock: sys.stderr.write("{}".format(msg))


def logArgs(args, width):
    text = "\n"
    for key in args:
        text += "\t{}: {}\n".format(key, args[key])
    log(msg=text, level="INFO")


def error(msg):
    return "{}{}{}".format(bcolors.FAIL, msg, bcolors.ENDC)


def close(msg):
    log(msg=msg, level="WARN")
    sys.exit(0)


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    BBLUE = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
