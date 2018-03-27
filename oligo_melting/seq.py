# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for basic nucleic acid sequence management.
'''

# CONSTANTS ====================================================================

AB_DNA = ["ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN"]
AB_RNA = ["ACGURYKMSWBDHVN", "UGCAYRMKSWVHDBN"]

# FUNCTIONS ====================================================================

def rc(na, t):
    '''
    Args:
        na (string): nucleic acid sequence.
        t (string): nucleic acid type, either 'dna' or 'rna'.

    Return:
        string: reverse complement of na.
    '''

    # Identify type
    t = t.lower()

    # Select alphabet
    if t == 'dna':
        ab = AB_DNA
    elif t == 'rna':
        ab = AB_RNA
    else:
        print('ERROR: unknown na type.')
        return()

    rab = ab[1].strip().lower()
    ab = ab[0].strip().lower()

    # Check provided string
    na = na.lower()

    for c in na:
        if not c in ab:
            print('ERROR: provided string conflicts with selected alphabet.')
            return()

    # Calculate reverse
    r = na[::-1]

    # Calculate reverse complement
    rc = []
    for c in r:
        rc.append(rab[ab.index(c)])

    rc=''.join([str(c) for c in rc]).upper()

    return(rc)

# END ==========================================================================

################################################################################
