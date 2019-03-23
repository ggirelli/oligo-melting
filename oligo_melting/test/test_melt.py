
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
'''

import oligo_melting as om


def test_seq():
    assert om.Sequence.check_ab("ACTG", om.AB_NA["DNA"])
    assert om.Sequence.check_ab("ACUG", om.AB_NA["RNA"])
    assert "CAGT" == om.Sequence.rc("ACTG", "DNA")
    assert "CUGT" == om.Sequence.rc("ACUG", "RNA")
    assert ["CU", "UG", "GT"] == om.Sequence.dimerator("CUGT")

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
