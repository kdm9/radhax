from __future__ import print_function, division, absolute_import
from collections import Counter
from argparse import ArgumentParser, FileType
import sys

import numpy as np
import screed

from .digest import Digest
from .utils import (
    clamp,
    output_frag_fasta,
    output_bed,
    perror,
    seqfile_iter_frags,
)


def add_common_args(ap):
    ap.add_argument('--genome', '-i', required=True, type=str,
                    help='Genome sequence (in fasta format)')
    ap.add_argument('--enzyme', '-e', required=True,
                    help='Restriction enzyme name')
    ap.add_argument('--enzyme2', '-r', default=None,
                    help='Second restriction enzyme name (for ddRAD, etc)')


def add_frag_len_args(ap):
    ap.add_argument('--max', '-x', default=10000, type=int,
                    help='Maximum fragment size')
    ap.add_argument('--min', '-n', default=1, type=int,
                    help='Minimum fragment size')


def hist_main():
    ap = ArgumentParser(description="Create histogram of fragment sizes")
    add_common_args(ap)
    add_frag_len_args(ap)
    ap.add_argument('--output', '-o', type=FileType('w'), default=sys.stdout,
                    help='Output file (default stdout)')
    ap.add_argument('--bins', '-b', type=int, default=100,
                    help='Number of bins in histogram')
    args = ap.parse_args()
    digestor = Digest(args.enzyme, args.enzyme2)

    sizes = []
    frags = seqfile_iter_frags(args.genome, digestor, minlen=args.min,
                               maxlen=args.max)
    for _, frag in frags:
        sizes.append(frag.len)
    sizes = np.array(sizes)
    counts, edges = np.histogram(sizes)
    last = edges[0]
    for i, count in enumerate(counts):
        range_str = "{}-{}".format(last, edges[i+1]-1)
        print(i+1, range_str, count, sep='\t', file=args.output)
        last = edges[i]


def digest_main():
    ap = ArgumentParser(description="Performs in-silico digestion of a genome")
    ap.add_argument('--output-fasta', type=FileType('w'),
                    help="Fasta file to output fragment sequences")
    ap.add_argument('--output-bed', type=FileType('w'),
                    help="Bed file to output fragment sequences")
    ap.add_argument('--ddrad', action="store_true",
                    help="Enforce different enzymes on each end of the frament.")
    add_common_args(ap)
    add_frag_len_args(ap)
    args = ap.parse_args()
    if not (args.output_fasta or args.output_bed):
        ap.error("One of --output-fasta FILE or --output-bed FILE is required")
    digestor = Digest(args.enzyme, args.enzyme2)

    frags = seqfile_iter_frags(args.genome, digestor, minlen=args.min,
                               maxlen=args.max, force_different_enzymes=args.ddrad)
    for read, frag in frags:
        if args.output_fasta:
            output_frag_fasta(read, frag, args.output_fasta)
        if args.output_bed:
            site_name = '{}[{}]_{}[{}]'.format(frag.lhs, frag.lhs_enzyme,
                                               frag.rhs, frag.rhs_enzyme)
            output_bed(read.name, frag.lhs, frag.rhs, site_name, args.output_bed)


def rebed_main():
    ap = ArgumentParser(description="Produces a BED file containing RE sites")
    ap.add_argument('--length', '-l', type=int, default=0, metavar='INT',
                    help='Include a window of INT bases around each RE site in the output')
    ap.add_argument('--output', '-o', type=FileType('w'), default=sys.stdout,
                    help='Output file (default stdout)')
    add_common_args(ap)
    args = ap.parse_args()
    digest = Digest(args.enzyme, args.enzyme2)

    for read in screed.open(args.genome, parse_description=True):
        seq = read.sequence
        seql = len(read.sequence)
        for cut, enz in digest.re_sites(seq):
            output_bed(read.name, clamp(cut-args.length, 0, seql),
                    clamp(cut + enz.size+args.length, 0, seql), str(enz), args.output)
