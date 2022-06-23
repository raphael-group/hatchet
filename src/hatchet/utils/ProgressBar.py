import sys
import datetime
import hatchet.utils.Supporting as sp


class ProgressBar:
    def __init__(
        self,
        total,
        length,
        counter=0,
        verbose=False,
        decimals=1,
        fill=chr(9608),
        lock=None,
        prefix='Progress:',
        suffix='Complete',
    ):
        self.total = total
        self.length = length
        self.decimals = decimals
        self.fill = fill
        self.prefix = prefix
        self.suffix = suffix
        self.lock = lock
        self.counter = counter
        self.verbose = verbose

    def progress(self, advance=True, msg=''):
        if self.lock is not None:
            self.progressLock(advance, msg)
        else:
            self.progressNoLock(advance, msg)

    def progressLock(self, advance=True, msg=''):
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            with self.counter.get_lock():
                self.counter.value += 1
        percent = ('{0:.' + str(self.decimals) + 'f}').format(100 * (self.counter.value / float(self.total)))
        filledLength = int(self.length * self.counter.value // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = '%s%s%s |%s| %s%s%% %s%s' % (
            sp.bcolors.BBLUE,
            self.prefix,
            sp.bcolors.ENDC,
            bar,
            sp.bcolors.BBLUE,
            percent,
            self.suffix,
            sp.bcolors.ENDC,
        )
        msg = '[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + msg
        if not self.verbose:
            toprint = rewind + result + ' %s[%s]%s' % (sp.bcolors.BBLUE, msg, sp.bcolors.ENDC)
        else:
            toprint = rewind + sp.bcolors.BBLUE + msg + sp.bcolors.ENDC + '\n' + result
        with self.lock:
            write(toprint)
            flush()
            if self.counter.value == self.total:
                write('\n')
                flush()

    def progressNoLock(self, advance=True, msg=''):
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            self.counter += 1
        percent = ('{0:.' + str(self.decimals) + 'f}').format(100 * (self.counter / float(self.total)))
        filledLength = int(self.length * self.counter // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = '%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix)
        if not self.verbose:
            toprint = rewind + result + ' [%s]' % (msg)
        else:
            toprint = rewind + msg + '\n' + result
        write(toprint)
        flush()
        if self.counter == self.total:
            write('\n')
            flush()
