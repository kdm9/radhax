from __future__ import print_function
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import Restriction
from Bio.Restriction import Analysis, RestrictionBatch
from collections import namedtuple
import screed
import sys

# Get all enzymes supported by the Bio.Retriction module
RE_ENZYMES = set(Restriction.Restriction_Dictionary.rest_dict.keys())

def list_enzymes(stream=sys.stderr):
    print("The following enzymes are supported:", file=stream)
    for enzyme_str in RE_ENZYMES:
        enzyme = getattr(Restriction, enzyme_str)
        print(enzyme, enzyme.site, sep='\t', file=stream)


# Container for fragment iteration
Fragment = namedtuple('Fragment', ['lhs', 'rhs', 'lhs_enzyme', 'rhs_enzyme'])


class Digest(object):
    '''Class whose methods digest sequences, returning different formats'''

    def __init__(self, enzyme, r2_enzyme=None):
        # If we don't have an r2 enzyme, use the r1 enzyme
        if r2_enzyme is None:
            r2_enzyme = enzyme
        self.enzyme = enzyme
        self.r2_enzyme = r2_enzyme
        self.enzyme_set = list(set([enzyme, r2_enzyme]))

    def iter_fragments(self, sequence, strict=True):
        '''Digests ``sequence``, and returns all fragments bordered by sites.

        if r2_enzyme is given, then only sites with ``enzyme`` at one end and
        ``r2_enzyme`` at the other are returned. Otherwise, all fragments are
        returned. If ``strict`` is False, all fragments are returned in any
        case.

        Returns fragments as a ``Fragment``.

        lhs/rhs are python slice intervals, i.e.:
        (first to include, first not to include)
        '''
        # TODO: enforce strictness
        seq = Seq(sequence, IUPACAmbiguousDNA)
        # Set up analysis class with our enzymes and seq
        rb = RestrictionBatch(self.enzyme_set)

        # Do digest and reformat to dict of {site: enz, site:enz}
        re_sites = {}
        for enzyme, cutsites in rb.search(seq).items():
            for cut in cutsites:
                re_sites[cut] = enzyme

        # loop through sites, yielding a Fragment for each
        last_enzyme = None
        this_pair = None
        for cut, enzyme in sorted(re_sites.items()):
            # Special case for the first site
            if this_pair is None:
                this_pair = [None, cut + enzyme.fst3 - 1 + enzyme.size]
                last_enzyme = enzyme
                continue
            # Switch time! this_pair[1] is the last cut site now.

            # The first base of the previous RE site (which was stored as the
            # rhs of a slice)
            this_pair[0] = this_pair[1] - last_enzyme.size
            # The first base not in the current RE site (rhs of the slice)
            this_pair[1] = cut + enzyme.fst3 - 1 + enzyme.size

            # Create fragment before switch below
            fragment = Fragment(lhs=this_pair[0], rhs=this_pair[1],
                                lhs_enzyme=last_enzyme, rhs_enzyme=enzyme)
            last_enzyme = enzyme
            yield fragment


def output_frag_fasta(read, frag, stream, width=80):
    print('>', read.name, '_', frag.lhs, '_', frag.rhs, sep='', file=stream)
    seq = read.sequence[frag.lhs:frag.rhs]
    for start in range(0, len(seq), width):
        print(seq[start:start+width], file=stream)


def output_frag_bed(read, frag, stream):
    site_name = '{}[{}]_{}[{}]'.format(frag.lhs, frag.lhs_enzyme, frag.rhs,
                                         frag.rhs_enzyme)
    print(read.name, frag.lhs, frag.rhs, site_name, sep='\t', file=stream)


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
