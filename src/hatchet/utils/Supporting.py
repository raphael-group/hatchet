import sys
import os
import os.path
import subprocess
import datetime
import requests
from urllib.parse import urlparse
import tarfile
from string import ascii_uppercase


def naturalOrder(text):
    return [int(s) if s.isdigit() else ord(s) for s in text]


def numericOrder(text):
    if len(digits(text)) > 0:
        return int(digits(text))
    else:
        if not (text.endswith('X') or text.endswith('Y')):
            raise ValueError(f'Found chromosome ID that is not numeric or X/Y: {text}')
        return ascii_uppercase.index(text[-1])


def digits(string):
    return "".join(x for x in string if x.isdigit())


def argmax(d):
    return max(d, key=(lambda x : d[x]))


def argmin(d):
    return min(d, key=(lambda x : d[x]))


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def url_exists(path):
    r = requests.head(path)
    return r.status_code == requests.codes.ok


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


# format strings with 'timestamp' and 'msg' as placeholders
MSG_FORMAT_STRINGS = {
    'STEP': bcolors.BOLD + bcolors.HEADER + '[{timestamp}]{msg}' + bcolors.ENDC,
    'INFO': bcolors.OKGREEN + '{msg}' + bcolors.ENDC,
    'WARN': bcolors.WARNING + '{msg}' + bcolors.ENDC,
    'PROGRESS': bcolors.BBLUE + '{msg}' + bcolors.ENDC,
    'ERROR': bcolors.FAIL + '{msg}' + bcolors.ENDC
}


def log(msg, level=None, lock=None, raise_exception=False, exception_class=ValueError):
    timestamp = '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
    format_string = MSG_FORMAT_STRINGS.get(level) or '{msg}'
    formatted_msg = format_string.format(msg=msg, timestamp=timestamp)

    if lock is None:
        sys.stderr.write(formatted_msg)
    else:
        with lock:
            sys.stderr.write(formatted_msg)

    if level == 'ERROR' and raise_exception:
        raise exception_class(msg)


def logArgs(args, width=40):
    text = '\n'
    for key in args:
        text += '\t{}: {}\n'.format(key, args[key])
    log(msg=text, level='INFO')


def error(msg, raise_exception=False, exception_class=ValueError):
    return log(msg, level='ERROR', raise_exception=raise_exception, exception_class=exception_class)


def ensure(pred, msg, exception_class=ValueError):
    if not pred:
        return error(msg, raise_exception=True, exception_class=exception_class)


def close(msg):
    log(msg=msg, level='WARN')
    sys.exit(0)


def run(commands, stdouterr_filepath=None, check_return_codes=True, error_msg=None, stdouterr_filepath_autoremove=True,
        **kwargs):
    singleton = isinstance(commands, str)
    if singleton:
        commands = [commands]
    return_codes = []

    # Very commonly used options in our code
    if 'universal_newlines' not in kwargs:
        kwargs['universal_newlines'] = True
    if 'shell' not in kwargs:
        kwargs['shell'] = True

    f = open(stdouterr_filepath, 'w') if stdouterr_filepath is not None else None
    for command in commands:
        p = subprocess.run(command, stdout=f, stderr=f, **kwargs)
        return_codes.append(p.returncode)
    if f:
        f.close()

    if check_return_codes:
        if error_msg is None:
            error_msg = 'Command "{command}" did not run successfully.'
        if stdouterr_filepath:
            error_msg += f' Please check {stdouterr_filepath} for possible hints.'

        for command, return_code in zip(commands, return_codes):
            if return_code != 0:
                error(error_msg.format(command=command), raise_exception=True)

    if stdouterr_filepath and stdouterr_filepath_autoremove:
        os.remove(stdouterr_filepath)

    return return_codes[0] if singleton else return_codes


def download(url, dirpath, overwrite=False, extract=True, sentinel_file=None, chunk_size=8192):
    out_basename = os.path.basename(urlparse(url).path)
    out_filepath = os.path.join(dirpath, out_basename)

    if sentinel_file is None or not os.path.exists(os.path.join(dirpath, sentinel_file)):
        if overwrite or not os.path.exists(out_filepath):
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(out_filepath, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        f.write(chunk)

        if (out_basename.endswith('.tar.gz') or out_basename.endswith('.tgz')) and extract:
            t = tarfile.open(out_filepath)
            t.extractall(dirpath)
            t.close()

    return out_filepath
