#!/usr/bin/env python
# encoding: utf-8


import sys
import errno
import argparse
import logging
from gzip import GzipFile
from bz2 import BZ2File
from collections import defaultdict, Counter

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
        Finds insertion points of mapped sequences.
        Expectation is that the sequences are upstream flanking the 5' of a insertion sequence.
        The query sequence is *always* as per the forward strand of the insert,
        to make deciding the start position of the insert easy.
        So, given a read, the following possibilities exist:
        1. If the read mapped to the forward strand of the insert, the flank is:
            >>>>>>>>>>>>>>
            [flank]^insert

        2. If the read mapped to the reverse of the insert, the flank is
            <<<<<<<<<<<<<<
            insert^[flank]
            and therefore, the extracted sequence is reverse-complemented

        This means that after mapping to the genome reference,
        1. If the flank maps to the forward strand, the insertion point is 3' (at the end) of the mapped sequence
        2. If it maps to the reverse, the insertion point is 5' of the mapped sequence (relative to the reference)
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
        # TODO add options here
        p.add_argument('files', metavar='file', nargs='*', help='SAM/BAM File to read. If none, read from STDIN')
        self.options = p.parse_args(args)

    def main(self):

        files = self.options.files
        if not files:
            files = '-'
        insertion_points = defaultdict(Counter)
        for f in files:
            with pysam.AlignmentFile(f, 'rb') as sam:
                for read in sam:
                    if not read.is_unmapped:
                        loc = None
                        if read.is_reverse:  # TODO: this might be off by one on either strand
                            loc = (read.reference_start, read.reference_start+1)
                        else:
                            loc = (read.reference_end, read.reference_end+1)

                        insertion_points[read.reference_id][loc] += 1

                print 'track type=bedGraph name="ura3 insertions"'
                for ref_id in sorted(insertion_points.keys()):
                    # with open('{}.bedGraph'.format(sam.references[ref_id]), 'w') as out
                    # print out 'browser position {}:0-{}'.format(sam.references[ref_id],sam.lengths[ref_id]-1)
                    for loc in sorted(insertion_points[ref_id].keys()):
                        count = insertion_points[ref_id][loc]
                        print '{}\t{}\t{}\t{}'.format(sam.references[ref_id], loc[0], loc[1], count)


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
