import radsim
import six


def test_list_enzymes():
    '''Check radsim.list_enzymes doesn't raise, and prints something'''
    sio = six.StringIO()
    radsim.list_enzymes(stream=sio)
    assert len(sio.getvalue()) > 0, sio.getvalue()


def test_digest_re_sites():
    '''Test the behaviour of Digest.re_sites()'''
    from radsim import Digest, Fragment
    from Bio.Restriction import PstI
    dig = Digest(PstI)
    seq = "NNCTGCAGacgtaCTGCAGacgtCTGCAGNN"
    #      0123456789012345678901234567890
    #        l          L     r         R
    # Above, lowercase is first frag, uppercase is 2nd frag
    got = list(dig.re_sites(seq))
    expected = [(2, PstI), (13, PstI), (23, PstI), ]
    assert got == expected


def test_digest_iter_fragments_simple():
    '''Test the behaviour of Digest.iter_fragments()'''
    from radsim import Digest, Fragment
    from Bio.Restriction import PstI
    pst1 = Digest(PstI)
    seq = "NNCTGCAGacgtaCTGCAGacgtCTGCAGNN"
    #      0123456789012345678901234567890
    #        l          L     r         R
    # Above, lowercase is first frag, uppercase is 2nd frag
    frags = list(pst1.iter_fragments(seq))
    expected = [Fragment(2, 19, PstI, PstI, 19 - 2),
                Fragment(13, 29, PstI, PstI, 29 - 13)]

    assert frags == expected, frags
    seq1 = seq[frags[0].lhs:frags[0].rhs]
    seq2 = seq[frags[1].lhs:frags[1].rhs]
    assert seq1 == 'CTGCAGacgtaCTGCAG', seq1
    assert seq2 == 'CTGCAGacgtCTGCAG', seq2


def test_digest_iter_fragments_2_enzymes():
    '''Test the behaviour of Digest.iter_fragments() w/ 2 enzymes'''
    from radsim import Digest, Fragment
    from Bio.Restriction import PstI, EcoRV
    pst1_ecor5 = Digest(PstI, EcoRV)
    seq = "NNCTGCAGacgtaGATATCacgtCTGCAGNN"
    #      0123456789012345678901234567890
    #        l          L     r         R
    # Above, lowercase is first frag, uppercase is 2nd frag
    frags = list(pst1_ecor5.iter_fragments(seq))
    expected = [Fragment(2, 19, PstI, EcoRV, 19 - 2),
                Fragment(13, 29, EcoRV, PstI, 29 - 13)]

    assert frags == expected, frags
    seq1 = seq[frags[0].lhs:frags[0].rhs]
    seq2 = seq[frags[1].lhs:frags[1].rhs]
    assert seq1 == 'CTGCAGacgtaGATATC', seq1
    assert seq2 == 'GATATCacgtCTGCAG', seq2
