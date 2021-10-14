from __future__ import print_function, division, absolute_import
from Bio.Seq import Seq
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
Fragment = namedtuple('Fragment', ['lhs', 'rhs', 'lhs_enzyme', 'rhs_enzyme',
                                   'len'])


class Digest(object):
    '''Class whose methods digest sequences, returning different formats'''

    def __init__(self, enzyme, r2_enzyme=None):
        # If we don't have an r2 enzyme, use the r1 enzyme
        if r2_enzyme is None:
            r2_enzyme = enzyme
        self.enzyme = enzyme
        self.r2_enzyme = r2_enzyme
        self.enzyme_set = list(set([enzyme, r2_enzyme]))

    def re_sites(self, sequence):
        seq = Seq(sequence)
        # Set up analysis class with our enzymes and seq
        rb = RestrictionBatch(self.enzyme_set)

        # Do digest and reformat to dict of {site: enz, site:enz}
        re_sites = {}
        for enzyme, cutsites in rb.search(seq).items():
            for cut in cutsites:
                cut = cut + enzyme.fst3 - 1
                re_sites[cut] = enzyme
        return sorted(re_sites.items())

    def iter_fragments(self, sequence, force_different_enzymes=True, minlen=0,
                       maxlen=sys.maxsize):
        '''Digests ``sequence``, and returns all fragments bordered by sites.

        if r2_enzyme is given, then only sites with ``enzyme`` at one end and
        ``r2_enzyme`` at the other are returned. Otherwise, all fragments are
        returned. If ``force_different_enzymes`` is False, all fragments are
        returned in any case.

        Returns fragments as a ``Fragment``.

        lhs/rhs are python slice intervals, i.e.:
        (first to include, first not to include)
        '''

        # loop through sites, yielding a Fragment for each
        last_enzyme = None
        last_start = None
        for cut, enzyme in self.re_sites(sequence):
            # Special case for the first site
            if last_enzyme is None:
                last_start = cut
                last_enzyme = enzyme
                continue

            this_end = cut + enzyme.size

            fraglen = this_end - last_start
            if fraglen < minlen or fraglen > maxlen:
                continue
            if force_different_enzymes and self.r2_enzyme is not None and self.enzyme != self.r2_enzyme and last_enzyme == enzyme:
                continue
            # Create fragment before switch below
            fragment = Fragment(lhs=last_start, rhs=this_end, len=fraglen,
                                lhs_enzyme=last_enzyme, rhs_enzyme=enzyme)
            last_enzyme = enzyme
            last_start = cut
            yield fragment
