# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for melting temperature calculation and correction for
              oligonucleotide duplexes.
'''

# DEPENDENCIES =================================================================

import math
import oligo_melting as OligoMelt
import os
import re
import sys

# CONSTANTS ====================================================================

# Gas constant
R = 1.987 / 1000    # kcal / (K mol)

# Thermodynamic tables ---------------------------------------------------------

# Table from Freier et al, PNAS(83), 1986 - in 1 M NaCl [RNA]
# dH0: kcal / mol
# dS0: eu = cal / (K mol)
# dG0: kcal / mol
NN_TABLE_RNA_RNA = {
    #                dH0     dS0     dG0
    'AA'        :   (-6.6,  -18.4,  -0.9),
    'UU'        :   (-6.6,  -18.4,  -0.9),
    'AU'        :   (-5.7,  -15.5,  -0.9),
    'UA'        :   (-8.1,  -22.6,  -1.1),
    'CA'        :   (-10.5, -27.8,  -1.8),
    'UG'        :   (-10.5, -27.8,  -1.8),
    'TG'        :   (-10.5, -27.8,  -1.8),
    'CU'        :   (-7.6,  -19.2,  -1.7),
    'AG'        :   (-7.6,  -19.2,  -1.7),
    'GA'        :   (-13.3, -35.5,  -2.3),
    'UC'        :   (-13.3, -35.5,  -2.3),
    'GU'        :   (-10.2, -26.2,  -2.1),
    'AC'        :   (-10.2, -26.2,  -2.1),
    'CG'        :   (-8.0,  -19.4,  -2.0),
    'GC'        :   (-14.2, -34.9,  -3.4),
    'GG'        :   (-12.2, -29.7,  -2.9),
    'CC'        :   (-12.2, -29.7,  -2.9),
    'endA'      :   (0,     -10.8,  3.4),
    'endC'      :   (0,     -10.8,  3.4),
    'endG'      :   (0,     -10.8,  3.4),
    'endU'      :   (0,     -10.8,  3.4),
    'has_end'   :   True,
    'has_init'  :   False,
    'sym'       :   (0,     -1.4,   0.4),
    'has_sym'   :   True
}

# Table from Sugimoto et al, Biochemistry(34), 1995 - in 1 M NaCl [DNA/RNA dupl]
# For calculation from the RNA sequence 5'-to-3'
# dH0: kcal / mol
# dS0: cal / (K mol)
# dG0: kcal / mol
NN_TABLE_RNA_DNA = {
    #               dH0     dS0     dG0
    'AA'        :   (-7.8,  -21.9,  -1.0),
    'AC'        :   (-5.9,  -12.3,  -2.1),
    'AG'        :   (-9.1,  -23.5,  -1.8),
    'AU'        :   (-8.3,  -23.9,  -0.9),
    'CA'        :   (-9.0,  -26.1,  -0.9),
    'CC'        :   (-9.3,  -23.2,  -2.1),
    'CG'        :   (-16.3, -47.1,  -1.7),
    'CU'        :   (-7.0,  -19.7,  -0.9),
    'GA'        :   (-5.5,  -13.5,  -1.3),
    'GC'        :   (-8.0,  -17.1,  -2.7),
    'GG'        :   (-12.8, -31.9,  -2.9),
    'GU'        :   (-7.8,  -21.6,  -1.1),
    'UA'        :   (-7.8,  -23.2,  -0.6),
    'UC'        :   (-8.6,  -22.9,  -1.5),
    'UG'        :   (-10.4, -28.4,  -1.6),
    'UU'        :   (-11.5, -36.4,  -0.2),
    'has_end'   :   False,
    'init'      :   (1.9,   -3.9,   3.1),
    'has_init'  :   True,
    'has_sym'   :   False
}

# Table from Sugimoto et al, Biochemistry(34), 1995 - in 1 M NaCl [DNA/RNA dupl]
# For calculation from the DNA sequence 5'-to-3'
# dH0: kcal / mol
# dS0: cal / (K mol)
# dG0: kcal / mol
NN_TABLE_DNA_RNA = {
    #               dH0     dS0     dG0
    'TT'        :   (-7.8,  -21.9,  -1.0),
    'GT'        :   (-5.9,  -12.3,  -2.1),
    'CT'        :   (-9.1,  -23.5,  -1.8),
    'AT'        :   (-8.3,  -23.9,  -0.9),
    'TG'        :   (-9.0,  -26.1,  -0.9),
    'GG'        :   (-9.3,  -23.2,  -2.1),
    'CG'        :   (-16.3, -47.1,  -1.7),
    'AG'        :   (-7.0,  -19.7,  -0.9),
    'TC'        :   (-5.5,  -13.5,  -1.3),
    'GC'        :   (-8.0,  -17.1,  -2.7),
    'CC'        :   (-12.8, -31.9,  -2.9),
    'AC'        :   (-7.8,  -21.6,  -1.1),
    'TA'        :   (-7.8,  -23.2,  -0.6),
    'GA'        :   (-8.6,  -22.9,  -1.5),
    'CA'        :   (-10.4, -28.4,  -1.6),
    'AA'        :   (-11.5, -36.4,  -0.2),
    'init'      :   (1.9,   -3.9,   3.1),
    'has_init'  :   True,
    'has_end'   :   False,
    'has_sym'   :   False
}

# Table from Allawi&Santalucia, Biochemistry(36), 1997 - in 1 M NaCl [DNA]
# dH0: kcal / mol
# dS0: eu = cal / (K mol)
# dG0: kcal / mol
NN_TABLE_DNA_DNA = {
    #                dH0     dS0     dG0
    'AA'        :   (-7.9,  -22.2,  -1.0),
    'TT'        :   (-7.9,  -22.2,  -1.0),
    'AT'        :   (-7.2,  -20.4,  -0.88),
    'TA'        :   (-7.2,  -21.3,  -0.58),
    'CA'        :   (-8.5,  -22.7,  -1.45),
    'TG'        :   (-8.5,  -22.7,  -1.45),
    'GT'        :   (-8.4,  -22.4,  -1.44),
    'AC'        :   (-8.4,  -22.4,  -1.44),
    'CT'        :   (-7.8,  -21.0,  -1.28),
    'AG'        :   (-7.8,  -21.0,  -1.28),
    'GA'        :   (-8.2,  -22.2,  -1.3),
    'TC'        :   (-8.2,  -22.2,  -1.3),
    'CG'        :   (-10.6, -27.2,  -2.17),
    'GC'        :   (-9.8,  -24.4,  -2.24),
    'GG'        :   (-8.0,  -19.9,  -1.84),
    'CC'        :   (-8.0,  -19.9,  -1.84),
    'endG'      :   (.1,    -2.8,   .98),
    'endC'      :   (.1,    -2.8,   .98),
    'endA'      :   (2.3,   4.1,    1.03),
    'endT'      :   (2.3,   4.1,    1.03),
    'has_end'   :   True,
    'has_init'  :   False,
    'sym'       :   (2.3,   4.1,    1.03),
    'has_sym'   :   True
}

# Tables and labels
NN_TABLES = {
    "DNA:DNA" : NN_TABLE_DNA_DNA, "RNA:RNA" : NN_TABLE_RNA_RNA, 
    "DNA:RNA" : NN_TABLE_DNA_RNA, "RNA:DNA" : NN_TABLE_RNA_DNA
}
NN_LABELS = list(NN_TABLES.keys())
NN_HYB_LABELS = ["DNA:RNA", "RNA:DNA"]
NN_DNA_TEMPLATE_LABELS = ["DNA:RNA", "DNA:DNA"]
NN_RNA_TEMPLATE_LABELS = ["RNA:RNA", "RNA:DNA"]

# Formamide correction ---------------------------------------------------------

FA_MODE_LABELS = ["wright", "mcconaughy"]
FA_MODE_WRIGHT2014 = FA_MODE_LABELS[0]
FA_MODE_MCCONA1969 = FA_MODE_LABELS[1]

# Default values ---------------------------------------------------------------

DEFAULT_NN_LABEL = NN_LABELS[0]
DEFAULT_OLIGO_CONC = 0.25e-6
DEFAULT_NA_CONC = 50e-3
DEFAULT_MG_CONC = 0
DEFAULT_FA_CONC = 0
DEFAULT_FA_MODE = FA_MODE_WRIGHT2014
DEFAULT_FA_MVALUE = "0.1734"
DEFAULT_T_CURVE_RANGE = 20
DEFAULT_T_CURVE_STEP = 0.5
DEFAULT_OUT_CURVE = False

# FUNCTIONS ====================================================================

def set_default(d, v = None):
    '''Set variable default value.

    Args:
        d: default value.
        v: current value.

    Args:
        d if v is None.
    '''
    if type(None) == type(v): return(d)
    return(v)

def get_dimers(seq):
    '''Extract NN dimers from sequence.

    Args:
        seq (str): nucleic acid sequence.

    Returns:
        list: NN dimers.
    '''
    return([seq[i:(i+2)] for i in range(len(seq) - 1)])

def parse_mval_string(fa_mval_s, fa_mode):
    '''Parse formamide m-value string.

    Args:
        fa_mval_s (str): formamide m-value string.
        fa_mode (str): formamide mode (see FA_MODE_LABELS).

    Returns:
        str, fun: formamide m-value string and function for parsing.
    '''

    # Evaluate formamide m-value
    parsed_mval = False
    mregex = ["(([+-]?[0-9\.]*)L([+-][0-9\.]*))", "([+-]?[0-9\.]*)"]

    # Search for xL+y m-value format
    msearch = re.search(mregex[0], fa_mval_s)
    if not type(None) is type(msearch):
        mgroups = msearch.groups()
        if fa_mval_s == mgroups[0]:
            mvalue = lambda x: float(mgroups[1]) * x + float(mgroups[2])
            parsed_mval = True

    # Search for x m-value format
    msearch = re.search(mregex[1], fa_mval_s)
    if not parsed_mval and not type(None) is type(msearch):
        mgroup = msearch.groups()[0]
        if fa_mval_s == mgroup:
            mvalue = lambda x: float(mgroup)
            parsed_mval = True

    # Trigger error otherwise
    if not parsed_mval:
        msg = "!!!ERROR! Unexpected formamide m-value format Check help page."
        sys.exit(msg)

    # Fix m-value label
    if not fa_mode in [FA_MODE_WRIGHT2014]:
        fa_mval_s = "-"

    return((fa_mval_s, mvalue))

def adj_fa(tm, h, s, seq, oligo_conc,
    fa_conc, fa_mode, mvalue, tt_mode = None, fa_conc_0 = None):
    '''Adjust melting temperature of a duplex based on formamide concentration.
    Based on Wright, Appl. env. microbiol.(80), 2014
    Or on McConaughy, Biochemistry(8), 1969
    
    Args:
      tm (float): melting temperature.
      h (float): standard enthalpy, in kcal / mol.
      s (float): standard enthropy, in kcal / (K mol).
      seq (string): oligonucleotide sequence.
      oligo_conc (float): oligonucleotide concentration in M.
      fa_conc (float): formamide concentration in %v,v.
      fa_mode (string): formamide correction lavel.
      mvalue (lambda): formamide m-value function.
      tt_mode (string): thermodynamic table label.
      fa_conc_0 (float): initial formamide concentration in %v,v.
    
    Returns:
      float: corrected melting temperature.
    '''
    
    # Default values -----------------------------------------------------------

    # Set default values
    tt_mode = set_default(DEFAULT_NN_LABEL, tt_mode)

    # Start --------------------------------------------------------------------


    # Default initial formamide concentration
    if type(None) == type(fa_conc_0):
        fa_conc_0 = 0

    # Delta formamide concentration
    dfac = fa_conc - fa_conc_0

    # If no change
    if 0 == dfac:
        return(tm)
    
    # Apply linear to melting temperature
    if FA_MODE_MCCONA1969 == fa_mode:
        tm -= 0.72 * dfac

    # Apply to free energy
    if FA_MODE_WRIGHT2014 == fa_mode:
        m = mvalue(len(seq))
        if tt_mode in NN_HYB_LABELS:
            tm = (h + m * dfac) / (R * math.log(oligo_conc / 4) + s)
        else:
            tm = (h + m * dfac) / (R * math.log(oligo_conc) + s)

    return(tm)

def adj_na(tm, na_conc, fgc, na_conc_0 = None):
    '''Adjust melting temperature of a duplexx based on sodium concentration.
    Based on Owczarzy et al, Biochemistry(43), 2004
    
    Args:
      tm (float): melting temperature at [Na] = 1 M.
      na_conc (float): monovalent species concentration in M.
      fgc (float): GC content fraction.
      na_conc_0 (float): initial monovalent species concentration.
      
    Returns:
      float: adjusted melting temperature.
    '''
    
    # Default
    if type(None) is type(na_conc_0):
        na_conc_0 = 1.
    if na_conc == na_conc_0:
        return(tm)
    Tm2 = tm

    # Parameters from paper
    na_a = 4.29e-5
    na_b = 3.95e-5
    na_c = 9.4e-6

    # Adjust
    if 1 == na_conc:
        Tm2r = (1. / tm)
        Tm2r -= (na_a * fgc - na_b) * math.log(na_conc_0 / na_conc)
        Tm2r -= na_c * (math.log(na_conc_0 / na_conc) ** 2)
        Tm2 = 1. / Tm2r
    elif 0 < na_conc:
        Tm2r = (1. / tm)
        Tm2r += (na_a * fgc - na_b) * math.log(na_conc / na_conc_0)
        Tm2r += na_c * (math.log(na_conc / na_conc_0) ** 2)
        Tm2 = 1. / Tm2r

    return(Tm2)

def adj_mg(tm, mg_conc, seq):
    '''Adjust melting temperature a duplexx based on sodium concentration.
    Based on Owczarzy et al, Biochemistry(47), 2008
    
    Args:
      tm (float): melting temperature at [Mg] = 0 M.
      mg_conc (float): divalent species concentration in M.
      seq (str): input sequence.
      
    Returns:
      float: adjusted melting temperature.
    '''
    
    # Default
    Tm3 = tm

    # Parameters from paper
    mg_a = 3.92e-5
    mg_b = -9.11e-6
    mg_c = 6.26e-5 
    mg_d = 1.42e-5
    mg_e = -4.82e-4
    mg_f = 5.25e-4
    mg_g = 8.31e-5

    # Adjust
    fgc = (seq.count('G') + seq.count('C')) / float(len(seq))
    if 0 < mg_conc:
        mg_conc_log = math.log(mg_conc)
        Tm3r = 1./tm + mg_a + mg_b*mg_conc_log + fgc*(mg_c + mg_d*mg_conc_log)
        Tm3r_factor = mg_e + mg_f*mg_conc_log + mg_g*(mg_conc_log)**2
        Tm3r_factor *= (1./(2*(len(seq) - 1)))
        Tm3r += Tm3r_factor
        Tm3 = 1./Tm3r

    return(Tm3)

def adj_ions(tm, na_conc, mg_conc, seq, na_conc_0 = None):
    '''Adjust melting temperature a duplexx based on ion concentration
    
    Args:
      tm (float): melting temperature.
      na_conc (float): monovalent species concentration in M.
      mg_conc (float): divalent species concentration in M.
      seq (str): input sequence.
      
    Returns:
      float: adjusted melting temperature.
    '''

    fgc = (seq.count('G') + seq.count('C')) / float(len(seq))
    
    if type(None) == type(na_conc_0):
        na_conc_0 = 1.

    if 0 != mg_conc:
        return(adj_mg(tm, mg_conc, seq))
    elif 0 != na_conc:
        return(adj_na(tm, na_conc, fgc, na_conc_0))
    else:
        return(tm)

def melt_curve(seq, h, s, tm, fgc, mvalue, oligo_conc = None, na_conc = None,
    mg_conc = None, fa_conc = None, fa_mode = None, trange = None, tstep = None,
    tt_mode = None):
    '''Generate melting curve
    
    Args:
      seq (string): oligonucleotide sequence.
      oligo_conc (float): oligonucleotide concentration in M.
      na_conc (float): monovalent species concentration in M.
      mg_conc (float): divalent species concentration in M.
      fa_conc (float): formamide concentration in %v,v.
      fa_mode (string): formamide correction lavel.
      mvalue (lambda): formamide m-value function.
      fgc (float): GC content fraction.
      h (float): standard enthalpy, in kcal / mol.
      s (float): standard enthropy, in kcal / (K mol).
      tm (float): melting temperature.
      trange (float): melting curve temperature range.
      tstep (float): melting curve temperature step.
      tt_mode (string): thermodynamic table label.
      
    Returns:
      list: melting curve data (x, y)
    '''

    # Default values -----------------------------------------------------------

    # Set default values
    oligo_conc = set_default(DEFAULT_OLIGO_CONC, oligo_conc)
    na_conc = set_default(DEFAULT_NA_CONC, na_conc)
    mg_conc = set_default(DEFAULT_MG_CONC, mg_conc)
    fa_conc = set_default(DEFAULT_FA_CONC, fa_conc)
    fa_mode = set_default(DEFAULT_FA_MODE, fa_mode)
    trange = set_default(DEFAULT_T_CURVE_STEP, curve_range)
    tstep = set_default(DEFAULT_T_CURVE_RANGE, curve_step)
    tt_mode = set_default(DEFAULT_NN_LABEL, tt_mode)

    # Start --------------------------------------------------------------------
    
    # Empty list for melting table
    data = []

    # Adjust melting temperature
    if FA_MODE_WRIGHT2014 == fa_mode:
        tm = adj_fa(tm, h, s, seq, oligo_conc,
        	fa_conc, fa_mode, mvalue, tt_mode)

    # Explore the temperature range
    t = tm - trange / 2.
    while t <= tm + trange / 2.:
        if FA_MODE_MCCONA1969 == fa_mode:
            # Calculate dissociated fraction
            k = dissoc_fraction(h, t, s, 0, oligo_conc, tt_mode)

            # Adjust output temperature
            t_out = adj_fa(t, h, s, seq, oligo_conc,
            	fa_conc, fa_mode, mvalue, tt_mode)
            t_out = adj_ions(t_out, na_conc, mg_conc, seq)

        if FA_MODE_WRIGHT2014 == fa_mode:
            # Calculate current FA m-value
            m = mvalue(len(seq))

            # Calculate dissociated fraction
            k = dissoc_fraction(h, t, s, m * fa_conc, oligo_conc, tt_mode)

            # Adjust output temperature
            t_out = adj_ions(t, na_conc, mg_conc, seq)

        # Append melting data
        data.append((t_out, k))
        
        # Go to next temperature
        t += tstep

    # Return melting table
    return(data)

def dissoc_fraction(h, t, s, gplus, oligo_conc = None, tt_mode = None):
    '''Calculate duplex dissociated fraction at given temperature.
    
    Args:
      h (float): enthalpy.
      t (float): current temperature.
      s (float): enthropy.
      oligo_conc (float): oligo concentration in M.
      tt_mode (string): N-N thermodynamic approach.
      gplus (float): additional free energy by solvent denaturant.
    
    Returns:
      float: dissociated fraction.
    '''

    # Default values -----------------------------------------------------------

    # Set default values
    oligo_conc = set_default(DEFAULT_OLIGO_CONC, oligo_conc)
    tt_mode = set_default(DEFAULT_NN_LABEL, tt_mode)

    # Start --------------------------------------------------------------------
    
    # Calculate free energy
    dg = h - t * s + gplus

    # Calculate factor
    if tt_mode in NN_HYB_LABELS:
        factor = math.exp(-dg / (R * t)) * (oligo_conc / 4)
    else:
        factor = math.exp(-dg / (R * t)) * oligo_conc

    # Calculate fraction
    k = 1 - factor / (1 + factor)

    # Output
    return(k)

def calc_tm_std(seq, tt_mode = None, couples = None, oligo_conc = None):
    '''Calculate melting temperature of a duplex at standard 1 M NaCl
    (monovalent ions conc). Based on SantaLucia, PNAS(95), 1998
    
    Args:
      seq (string): oligonucleotide sequence.
      tt_mode (string): thermodynamic table label.
      couples (list): list of base dimers.
      oligo_conc (float): oligonucleotide concentration in M.
    
    Returns:
      tuple: hybridization enthalpy, enthropy and melting temperature.
    '''
    
    # Default values -----------------------------------------------------------

    # Set default values
    oligo_conc = set_default(DEFAULT_OLIGO_CONC, oligo_conc)
    tt_mode = set_default(DEFAULT_NN_LABEL, tt_mode)
    if type(None) == type(couples):
        couples = get_dimers(seq)

    # Build upon input ---------------------------------------------------------

    # Select thermodynamic table
    tt = NN_TABLES[tt_mode]

    # Start --------------------------------------------------------------------
    
    # Standard 1 M NaCl
    na_conc_0 = 1.

    # Calculate dH0(N-N)
    h = sum([tt[c][0] for c in couples])

    # Add initiation enthalpy
    if tt['has_end']:
        h += tt['end%s' % (seq[0],)][0]
        h += tt['end%s' % (seq[-1],)][0]
    if tt['has_init']:
        h += tt['init'][0]

    # Correct enthalpy for symmetry
    if tt['has_sym'] and seq == OligoMelt.rc(seq, 'dna'):
        h += tt['sym'][0]

    # Calculate dS0(N-N) in kcal / (K mol)
    s = sum([tt[c][1] for c in couples])

    # Add initiation enthropy
    if tt['has_end']:
        s += tt['end%s' % (seq[0],)][1]
        s += tt['end%s' % (seq[-1],)][1]
    if tt['has_init']:
        s += tt['init'][1]

    # Correct enthalpy for symmetry
    if tt['has_sym'] and seq == OligoMelt.rc(seq, 'dna'):
        s += tt['sym'][1]
    s /= 1e3

    # Calculate melting temperature in Celsius
    if tt_mode in NN_HYB_LABELS:
        tm = h / (s + R * math.log(oligo_conc / 4))
    else:
        tm = h / (s + R * math.log(oligo_conc))

    # Output
    return((h, s, tm))

def calc_tm(seq, name = None, oligo_conc = None, na_conc = None, mg_conc = None,
	fa_conc = None, fa_mode = None, fa_mval_s = None,
    tt_mode = None, celsius = None, is_verbose = None,
    curve_step = None, curve_range = None, curve_outpath = None, silent = None):
    '''Calculate melting temperature of provided oligo duplex sequence.
    
    Args:
      seq (string): oligonucleotide sequence.
      oligo_conc (float): oligonucleotide concentration in M.
      na_conc (float): monovalent species concentration in M.
      mg_conc (float): divalent species concentration in M.
      fa_conc (float): formamide concentration in %v,v.
      fa_mode (string): formamide correction lavel.
      fa_mval_s (string): formamide m-value string.
      tt_mode (string): thermodynamic table label.
      celsius (bool): convert K to degC.
      is_verbose (bool): be verbose.
      curve_step (float): melting curve temperature step.
      curve_range (float): melting curve temperature range.
      curve_outpath (string): melting curve output path.

    Returns:
      (name, g, h, s, tm, seq)
      name (str): sequence name.
      g (float): delta Gibbs free energy.
      h (float): delta enthalpy.
      s (float): delta entropy.
      tm (float): melting temperature after corrections.
      seq (str): input sequence.
    '''

    # Default values -----------------------------------------------------------

    # Set default values
    name = set_default("seq", name)
    oligo_conc = set_default(DEFAULT_OLIGO_CONC, oligo_conc)
    na_conc = set_default(DEFAULT_NA_CONC, na_conc)
    mg_conc = set_default(DEFAULT_MG_CONC, mg_conc)
    fa_conc = set_default(DEFAULT_FA_CONC, fa_conc)
    fa_mode = set_default(DEFAULT_FA_MODE, fa_mode)
    fa_mval_s = set_default(DEFAULT_FA_MVALUE, fa_mval_s)
    tt_mode = set_default(DEFAULT_NN_LABEL, tt_mode)
    celsius = set_default(False, celsius)
    is_verbose = set_default(False, is_verbose)
    curve_step = set_default(DEFAULT_T_CURVE_RANGE, curve_step)
    curve_range = set_default(DEFAULT_T_CURVE_STEP, curve_range)
    curve_outpath = set_default(DEFAULT_OUT_CURVE, curve_outpath)
    silent = set_default(True, silent)

    # Build upon input ---------------------------------------------------------

    # Remove output curve file if it exists
    do_curve = False != curve_outpath
    if do_curve and os.path.isfile(curve_outpath):
        os.remove(curve_outpath)

    # Parse formamide mvalue string
    fa_mval_s, mvalue = parse_mval_string(fa_mval_s, fa_mode)

    # Select thermodynamic table
    tt = NN_TABLES[tt_mode]

    # Start --------------------------------------------------------------------

    # Make string uppercase
    seq = seq.upper()

    # Check that all characters code for bases
    if not 0 == sum([1 for b in seq if b not in 'TCGAU']):
        print('The provided sequence contains non-nucleic acid characters.')
        return

    # Check that the correct nucleic acid type was provided
    if tt_mode in NN_DNA_TEMPLATE_LABELS:
        if 'U' in seq:
            print('The option -t %s requires a DNA sequence.' % (tt_mode,))
            return
    elif tt_mode in NN_RNA_TEMPLATE_LABELS:
        if 'T' in seq:
            print('The option -t %s requires a RNA sequence.' % (tt_mode,))
            return

    # Make NN couples
    couples = get_dimers(seq)

    # Calculate GC content
    fgc = (seq.count('G') + seq.count('C')) / float(len(seq))

    # 1 M NaCl case
    # Based on SantaLucia, PNAS(95), 1998
    # -----------------------------------
    (h, s, Tm1) = calc_tm_std(seq, tt_mode, couples, oligo_conc)

    # Adjust for FA
    # Based on Wright, Appl. env. microbiol.(80), 2014
    # Or on McConaughy, Biochemistry(8), 1969
    # ------------------------------------------------
    Tm2 = adj_fa(Tm1, h, s, seq, oligo_conc,
    	fa_conc, fa_mode, mvalue, tt_mode)

    # Adjusted for [Na]
    # Based on Owczarzy et al, Biochemistry(43), 2004
    # -----------------------------------------------
    Tm3 = adj_na(Tm2, na_conc, fgc)

    # Adjusted for Mg
    # Based on Owczarzy et al, Biochemistry(47), 2008
    # -----------------------------------------------
    if 0 < mg_conc:
        Tm4 = adj_mg(Tm2, mg_conc, seq)
    else:
        Tm4 = Tm3

    # Generate melting curves
    # -----------------------
    
    if do_curve:
        if 0 == len(os.path.basename(curve_outpath)):
            msg = "\n!!!ERROR! --out-curve should be a file, not a folder."
            msg += "\n          Path: %s" % curve_outpath
            sys.exit(msg)
        fout = open(curve_outpath, 'a+')
        tab = melt_curve(seq, h, s, Tm1, fgc, mvalue, oligo_conc, na_conc,
            mg_conc, fa_conc, fa_mode, curve_range, curve_step, tt_mode)
        for (t, k) in tab:
            if celsius:
                fout.write("%s\t%f\t%f\n" % (name, t - 273.15, k))
            else:
                fout.write("%s\t%f\t%f\n" % (name, t, k))
        fout.close()

    # Log output
    # ----------

    g = h - (37 + 273.15) * s

    if not is_verbose:
        if celsius:
            Tm4 -= 273.15
        output = (name, g, h, s, Tm4, seq)
        if not silent: print("%s\t%f\t%f\t%f\t%f\t%s" % output)
        return(output)
    else:
        output = """
             Oligo label : %s
          Oligo sequence : %s
              GC-content : %.0f%%
                 [oligo] : %.9f M
                   [Na+] : %f M
                  [Mg2+] : %f M
                    [FA] : %.1f%%
           FA correction : %s
              FA m-value : %s
             Duplex type : %s

             Melt. curve : %s
             Curve range : %f
              Curve step : %f
                  Output : %s

                     dH0 : %f kcal/mol 
                     dS0 : %f kcal/(K·mol)
                     dG0 : %f kcal/mol

             [Na+] = 1 M : Tm = %f K (= %f degC)

      [FA] = %.1f %%(v,v) : Tm = %f K (= %f degC)

      [Na+] = %.6f M : Tm = %f K (= %f degC)

     [Mg2+] = %.6f M : Tm = %f K (= %f degC)
     Note: Mg2+ correction overwrites Na+ correction.
        """ % (
            name, seq, fgc * 100,
            oligo_conc, na_conc, mg_conc,
            fa_conc, fa_mode, fa_mval_s, tt_mode,
            do_curve, curve_range, curve_step, curve_outpath,
            h, s, g,
            Tm1, Tm1 - 273.15,
            fa_conc, Tm2, Tm2 - 273.15,
            na_conc, Tm3, Tm3 - 273.15,
            mg_conc, Tm4, Tm4 - 273.15,
        )
        if not silent: print(output)
        return(output)

# END ==========================================================================

################################################################################
