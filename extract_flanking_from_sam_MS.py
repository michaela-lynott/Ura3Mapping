#!/usr/bin/env python
# encoding: utf-8
#written by Alex Richter and edited by Michaela Lynott for mapping URA3 insertion points along K.max genome

import sys
import errno
import argparse
import logging
from gzip import GzipFile
from bz2 import BZ2File
import pysam


log = logging.getLogger(__name__)


# utility functions
def open_any(filename):
    """Opens uncompressed or compressed files, returning the file object,
    or, if passed a filename of '-', returns STDIN"""
    if filename == '-':
        fh = sys.stdin
    elif filename[-3:] == '.gz':
        fh = GzipFile(filename, 'r')
    elif filename[-4:] == '.bz2':
        fh = BZ2File(filename, 'r')
    else:
        fh = open(filename, 'r')

    return fh


class App(object):
    """
        Main script wrapper: handles argument parsing and file opening
    """

    def __init__(self, argv=None):
        if not argv:
            argv = sys.argv
        # parse arguments
        self.exe_name = argv[0]
        self.options = None
        self.parse_options(argv[1:])

        # set up the logger
        log_level = logging.WARNING
        if self.options.verbose:
            log_level = logging.DEBUG
        log.level = log_level

        log.debug('%s %s', self.exe_name, ' '.join(argv[1:]))

    def parse_options(self, args):
        p = argparse.ArgumentParser(description='<description of the program function>')
        p.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose logging')
        p.add_argument('-o', '--overlap', default=0, type=int, help='Overlap into mapped sequence')
        # TODO add options here
        p.add_argument('files', metavar='file', nargs='*', help='File to read. If none, read from STDIN')
        self.options = p.parse_args(args)

    def main(self):

        files = self.options.files
        if not files:
            files = '-'
        n = 0
        for f in files:
            with pysam.AlignmentFile(f, 'rb') as sam:
                for read in sam:  # only want 5' flanking region from insert perspective
                    n += 1
                    if (not read.is_unmapped) \
                            and (not read.is_secondary) \
                            and (not read.is_supplementary) \
                            and read.cigartuples[0][0] in (4, 5) \
                            and read.cigartuples[0][1] >= 30 \
                            and read.cigartuples[1][0] == 0 \
                            and read.reference_start <= 3:
                        unmapped_bases = read.cigartuples[0][1]
                        flank_seq = read.query_sequence[0:unmapped_bases+self.options.overlap]  # in sam file, seq is always ref-centric
                        if read.is_reverse:
                            q_end = read.query_length
                            q_start = q_end - unmapped_bases
                            print '>{} {} {}-{} revcomp\n{}\n'.format(n, read.query_name, q_start, q_end, flank_seq)
                        else:
                            print '>{} {} 0-{}\n{}\n'.format(n, read.query_name, unmapped_bases, flank_seq)
                    elif (not read.is_unmapped) \
                            and (not read.is_secondary) \
                            and (not read.is_supplementary) \
                            and read.cigartuples[-1][0] in (4,5) \
                            and read.cigartuples[-1][1] >= 30 \
                            and read.cigartuples[-2][0] == 0 \
                            and read.cigartuples[-2][1] >= 30:
                        unmapped_bases = read.cigartuples[-1][1]
                        flank_seq = read.query_sequence[-1*(unmapped_bases+self.options.overlap):]
                        if read.is_reverse:
                            q_end = read.query_length
                            q_start = q_end - unmapped_bases
                            print '>{} {} {}-{} revcomp\n{}\n'.format(n, read.query_name, q_start, q_end, flank_seq)
                        else:
                            print '>{} {} 0-{}\n{}\n'.format(n, read.query_name, unmapped_bases, flank_seq)
# -----------------------------------------
if __name__ == '__main__':
    logging.basicConfig()
    app = App()
    # trap Ctrl-C and sigpipe to block tracebacks
    try:
        retval = app.main()
    except KeyboardInterrupt:
        sys.exit(1)
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(1)
        else:
            raise
    else:
        sys.exit(retval)
