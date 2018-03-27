# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for melting temperature calculation and correction for
              single-strand oligonucleotide secondary structures.
'''

# DEPENDENCIES =================================================================

import math


# CONSTANTS ====================================================================

# Gas constant
R = 1.987 / 1000    # kcal / (K mol)

# FUNCTIONS ====================================================================

def adj_fa(tm, fa_conc, fa_conc_0 = None):
    '''Adjust secondary structure melting temperature based on formamide
    concentration. Method adapted to work with OligoArrayAux output.
    Based on McConaughy, Biochemistry(8), 1969
    Method based on Wright et al not implemented, as m-value is not available.
    
    Args:
      tm (float): melting temperature, in K.
      fa_conc (float): formamide concentration in %v,v.
      fa_mode (string): formamide correction lavel.
      fa_conc_0 (float): initial formamide concentration in %v,v.
    
    Returns:
      float: corrected melting temperature.
    '''
    
    # Default initial concentration
    if type(None) == type(fa_conc_0):
        fa_conc_0 = 0

    # Delta FA concentration
    dfac = fa_conc - fa_conc_0

    # Correct
    if 0 == dfac:
        return(tm)
    else:
        return(tm - 0.72 * dfac)

def melt_curve(h, s, tm, fa_conc, trange, tstep):
    '''Generate melting curve
    
    Args:
      tm (float): melting temperature.
      fa_conc (float): formamide concentration in %v,v.
      trange (float): melting curve temperature range.
      tstep (float): melting curve temperature step.
      
    Returns:
      list: melting curve data (x, y)
    '''

    # Empty list for melting table
    data = []

    # Explore the temperature range
    t = tm - trange / 2.
    while t <= tm + trange / 2.:
        # Calculate dissociated fraction
        k = unfolded_fraction(h, t, s)

        # Adjust output temperature
        t_out = adj_fa(t, fa_conc)

        # Append melting data
        data.append((t_out, k))

        # Go to next temperature
        t += tstep

    # Return melting table
    return(data)

def unfolded_fraction(h, t, s):
    '''Calculate secondary structure unfolded fraction
    
    Args:
      h (float): enthalpy.
      t (float): current temperature.
      s (float): enthropy.
    
    Returns:
      float: dissociated fraction.
    '''
    
    # Calculate free energy at current temperature
    dg = h - t * s

    # Calculate factor
    factor = math.exp(-dg / (R * t))

    # Calculate fraction
    k = 1 / (1 + factor)

    # Output
    return(k)

# END ==========================================================================

################################################################################
