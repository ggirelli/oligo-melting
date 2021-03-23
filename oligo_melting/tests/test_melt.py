"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from oligo_melting.const import AB_NA, NATYPES
from oligo_melting.sequence import Sequence


def test_seq():
    assert Sequence.check_ab("ACTG", AB_NA[NATYPES["DNA"]])
    assert Sequence.check_ab("ACUG", AB_NA[NATYPES["RNA"]])
    assert "CAGT" == Sequence("ACTG", NATYPES["DNA"]).rc
    assert "CAGU" == Sequence("ACUG", NATYPES["RNA"]).rc
    assert ["CU", "UG", "GT"] == [dimer for dimer in Sequence.dimerator("CUGT")]


def test_nnet():
    # Test NNET building class
    pass


def test_ions():
    # Test single monovalent correction
    # Test single divalent correction
    # Test workflow for ions correction
    pass


def test_denaturants():
    # Test McConaughy correction
    # Test mvalue calculation
    # Test Wright correction with m1 only
    # Test Wrihgt correction with m1 and m2
    pass


def test_melter():
    # Test standard ([Na+] = 1 M) Tm calculation
    # Test general calculation, with denaturant and mono/divalent ions
    # Test dissociated fraction calculation
    # Test melting curve generation
    pass
