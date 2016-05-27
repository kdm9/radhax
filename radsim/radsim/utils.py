from __future__ import print_function, division, absolute_import
import sys
import screed


def perror(*args, **kwargs):
    '''Print an error to stderr, and exit with non-zero code.'''
    print(file=sys.stderr, *args, **kwargs)
    exit(1)


def output_frag_fasta(read, frag, stream, width=80):
    print('>', read.name, '_', frag.lhs, '_', frag.rhs, sep='', file=stream)
    seq = read.sequence[frag.lhs:frag.rhs]
    for start in range(0, len(seq), width):
        print(seq[start:start+width], file=stream)


def output_frag_bed(read, frag, stream):
    site_name = '{}[{}]_{}[{}]'.format(frag.lhs, frag.lhs_enzyme, frag.rhs,
                                       frag.rhs_enzyme)
    print(read.name, frag.lhs, frag.rhs, site_name, sep='\t', file=stream)


def seqfile_iter_frags(seqfile, digestor, minlen, maxlen):
    for read in screed.open(seqfile):
        seq = read.sequence
        for frag in digestor.iter_fragments(seq, minlen=minlen, maxlen=maxlen):
            yield read, frag
