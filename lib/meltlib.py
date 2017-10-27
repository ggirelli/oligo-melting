# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.4.4
# Date: 20171026
# Project: oligo characterization
# Description: functions to calculate and correct melting temperature.
# 
# References:
#  [1] Freier et al, PNAS(83), 1986;
#  [2] Sugimoto et al, Biochemistry(34), 1995.
#  [3] Allawi & Santalucia, Biochemistry(36), 1997;
#  [4] SantaLucia, PNAS(95), 1998;
#  [5] Owczarzy et al, Biochemistry(43), 2004;
#  [6] Owczarzy et al, Biochemistry(47), 2008;
#  [7] McConaughy et al, Biochemistry(8), 1969;
#  [8] Wright et al, Appl. env. microbiol.(80), 2014.
# 
# ------------------------------------------------------------------------------

# DEPENDENCIES =================================================================

import math

# PARAMETERS ===================================================================

# Constants --------------------------------------------------------------------

# Gas constant
R = 1.987 / 1000    # kcal / (K mol)

# Thermodynamic tables ---------------------------------------------------------

# Table from Freier et al, PNAS(83), 1986 - in 1 M NaCl [RNA]
# dH0: kcal / mol
# dS0: eu = cal / (K mol)
# dG0: kcal / mol
freier = {
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
sugimotor = {
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
sugimotod = {
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
allawi = {
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
        ab = ["ATCG", "TAGC"]
    elif t == 'rna':
        ab = ["AUCG", "UAGC"]
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

def ssMelt_fa_adj(tm, fa_conc, fa_conc_0 = None):
    # Adjust secondary structure melting temperature based on formamide
    # concentration. Method adapted to work with OligoArrayAux output.
    # Based on McConaughy, Biochemistry(8), 1969
    # Method based on Wright et al not implemented, as m-value is not available.
    # 
    # Args:
    #   tm (float): melting temperature, in K.
    #   fa_conc (float): formamide concentration in %v,v.
    #   fa_mode (string): formamide correction lavel.
    #   fa_conc_0 (float): initial formamide concentration in %v,v.
    # 
    # Returns:
    #   float: corrected melting temperature.
    
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

def duMelt_fa_adj(tm, h, s, seq, oligo_conc,
    fa_conc, fa_mode, mvalue, tt_mode, fa_conc_0 = None):
    # Adjust melting temperature of a duplex based on formamide concentration.
    # Based on Wright, Appl. env. microbiol.(80), 2014
    # Or on McConaughy, Biochemistry(8), 1969
    # 
    # Args:
    #   tm (float): melting temperature.
    #   h (float): standard enthalpy, in kcal / mol.
    #   s (float): standard enthropy, in kcal / (K mol).
    #   seq (string): oligonucleotide sequence.
    #   oligo_conc (float): oligonucleotide concentration in M.
    #   fa_conc (float): formamide concentration in %v,v.
    #   fa_mode (string): formamide correction lavel.
    #   mvalue (lambda): formamide m-value function.
    #   tt_mode (string): thermodynamic table label.
    #   fa_conc_0 (float): initial formamide concentration in %v,v.
    # 
    # Returns:
    #   float: corrected melting temperature.
    
    # Default initial formamide concentration
    if type(None) == type(fa_conc_0):
        fa_conc_0 = 0

    # Delta formamide concentration
    dfac = fa_conc - fa_conc_0

    # If no change
    if 0 == dfac:
        return(tm)
    
    # Apply linear to melting temperature
    if "mcconaughy" == fa_mode:
        tm -= 0.72 * dfac

    # Apply to free energy
    if "wright" == fa_mode:
        m = mvalue(len(seq))
        if "sugimoto" in tt_mode:
            tm = (h + m * dfac) / (R * math.log(oligo_conc / 4) + s)
        else:
            tm = (h + m * dfac) / (R * math.log(oligo_conc) + s)

    return(tm)

def duMelt_na_adj(tm, na_conc, fgc, na_conc_0 = None):
    # Adjust melting temperature of a duplexx based on sodium concentration.
    # Based on Owczarzy et al, Biochemistry(43), 2004
    # 
    # Args:
    #   tm (float): melting temperature at [Na] = 1 M.
    #   na_conc (float): monovalent species concentration in M.
    #   fgc (float): GC content fraction.
    #   
    # Returns:
    #   float: adjusted melting temperature.
    
    # Default
    if type(None) is type(na_conc_0):
        na_conc_0 = 1.
    Tm2 = tm

    # Parameters from paper
    na_a = 4.29e-5
    na_b = 3.95e-5
    na_c = 9.4e-6

    # Adjust
    if 0 < na_conc:
        Tm2r = (1. / tm)
        Tm2r += (na_a * fgc - na_b) * math.log(na_conc / na_conc_0)
        Tm2r += na_c * (math.log(na_conc / na_conc_0) ** 2)
        Tm2 = 1. / Tm2r

    return(Tm2)

def duMelt_mg_adj(tm, mg_conc, fgc):
    # Adjust melting temperature a duplexx based on sodium concentration.
    # Based on Owczarzy et al, Biochemistry(47), 2008
    # 
    # Args:
    #   tm (float): melting temperature at [Mg] = 0 M.
    #   mg_conc (float): divalent species concentration in M.
    #   fgc (float): GC content fraction.
    #   
    # Returns:
    #   float: adjusted melting temperature.
    
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
    if 0 < mg_conc:
        mg_conc_log = math.log(mg_conc)
        Tm3r = 1./tm + mg_a + mg_b*mg_conc_log + fgc*(mg_c + mg_d*mg_conc_log)
        Tm3r_factor = mg_e + mg_f*mg_conc_log + mg_g*(mg_conc_log)**2
        Tm3r_factor *= (1./(2*(len(seq) - 1)))
        Tm3r += Tm3r_factor
        Tm3 = 1./Tm3r

    return(Tm3)

def duMelt_ion_adj(tm, na_conc, mg_conc, fgc):
    # Adjust melting temperature a duplexx based on ion concentration
    # 
    # Args:
    #   tm (float): melting temperature.
    #   na_conc (float): monovalent species concentration in M.
    #   mg_conc (float): divalent species concentration in M.
    #   fgc (float): GC content fraction.
    #   
    # Returns:
    #   float: adjusted melting temperature.
    if 0 != mg_conc:
        return(duMelt_mg_adj(tm, mg_conc, fgc))
    elif 0 != na_conc:
        return(duMelt_na_adj(tm, na_conc, fgc))
    else:
        return(tm)

def duMelt_curve(seq, oligo_conc, na_conc, mg_conc, fa_conc, fa_mode, mvalue,
    fgc, h, s, tm, trange, tstep, tt_mode):
    # Generate melting curve
    # 
    # Args:
    #   seq (string): oligonucleotide sequence.
    #   oligo_conc (float): oligonucleotide concentration in M.
    #   na_conc (float): monovalent species concentration in M.
    #   mg_conc (float): divalent species concentration in M.
    #   fa_conc (float): formamide concentration in %v,v.
    #   fa_mode (string): formamide correction lavel.
    #   mvalue (lambda): formamide m-value function.
    #   fgc (float): GC content fraction.
    #   h (float): standard enthalpy, in kcal / mol.
    #   s (float): standard enthropy, in kcal / (K mol).
    #   tm (float): melting temperature.
    #   trange (float): melting curve temperature range.
    #   tstep (float): melting curve temperature step.
    #   tt_mode (string): thermodynamic table label.
    #   
    # Returns:
    #   list: melting curve data (x, y)
    
    # Empty list for melting table
    data = []

    # Adjust melting temperature
    if "wright" == fa_mode:
        tm = duMelt_fa_adj(tm, h, s, seq, oligo_conc,
        	fa_conc, fa_mode, mvalue, tt_mode)

    # Explore the temperature range
    t = tm - trange / 2.
    while t <= tm + trange / 2.:
        if "mcconaughy" == fa_mode:
            # Calculate dissociated fraction
            k = duDissoc_fraction(h, t, s, oligo_conc, tt_mode, 0)

            # Adjust output temperature
            t_out = duMelt_ion_adj(t, na_conc, mg_conc, fgc)
            t_out = duMelt_fa_adj(t_out, h, s, seq, oligo_conc,
            	fa_conc, fa_mode, mvalue, tt_mode)

        if "wright" == fa_mode:
            # Calculate current FA m-value
            m = mvalue(len(seq))

            # Calculate dissociated fraction
            k = duDissoc_fraction(h, t, s, oligo_conc, tt_mode, m * fa_conc)

            # Adjust output temperature
            t_out = duMelt_ion_adj(t, na_conc, mg_conc, fgc)

        # Append melting data
        data.append((t_out, k))
        
        # Go to next temperature
        t += tstep

    # Return melting table
    return(data)

def ssMelt_curve(h, s, tm, fa_conc, trange, tstep):
    # Generate melting curve
    # 
    # Args:
    #   tm (float): melting temperature.
    #   fa_conc (float): formamide concentration in %v,v.
    #   trange (float): melting curve temperature range.
    #   tstep (float): melting curve temperature step.
    #   
    # Returns:
    #   list: melting curve data (x, y)

    # Empty list for melting table
    data = []

    # Explore the temperature range
    t = tm - trange / 2.
    while t <= tm + trange / 2.:
        # Calculate dissociated fraction
        k = ssUnfolded_fraction(h, t, s)

        # Adjust output temperature
        t_out = ssMelt_fa_adj(t, fa_conc)

        # Append melting data
        data.append((t_out, k))

        # Go to next temperature
        t += tstep

    # Return melting table
    return(data)

def duDissoc_fraction(h, t, s, oligo_conc, tt_mode, gplus):
    # Calculate duplex dissociated fraction at given temperature.
    # 
    # Args:
    #   h (float): enthalpy.
    #   t (float): current temperature.
    #   s (float): enthropy.
    #   oligo_conc (float): oligo concentration in M.
    #   tt_mode (string): N-N thermodynamic approach.
    #   gplus (float): additional free energy by solvent denaturant.
    # 
    # Returns:
    #   float: dissociated fraction.
    
    # Calculate free energy
    dg = h - t * s + gplus

    # Calculate factor
    if "sugimoto" in tt_mode:
        factor = math.exp(-dg / (R * t)) * (oligo_conc / 4)
    else:
        factor = math.exp(-dg / (R * t)) * oligo_conc

    # Calculate fraction
    k = 1 - factor / (1 + factor)

    # Output
    return(k)

def ssUnfolded_fraction(h, t, s):
    # Calculate secondary structure unfolded fraction
    # 
    # Args:
    #   h (float): enthalpy.
    #   t (float): current temperature.
    #   s (float): enthropy.
    # 
    # Returns:
    #   float: dissociated fraction.
    
    # Calculate free energy at current temperature
    dg = h - t * s

    # Calculate factor
    factor = math.exp(-dg / (R * t))

    # Calculate fraction
    k = 1 / (1 + factor)

    # Output
    return(k)

def duMelt_std_calc(seq, tt, tt_mode, couples, oligo_conc):
    # Calculate melting temperature of a duplex at standard 1 M NaCl (monovalent
    # ions conc). Based on SantaLucia, PNAS(95), 1998
    # 
    # Args:
    #   seq (string): oligonucleotide sequence.
    #   tt (dict): thermodynamic table list.
    #   tt_mode (string): thermodynamic table label.
    #   couples (list): list of base dimers.
    #   oligo_conc (float): oligonucleotide concentration in M.
    # 
    # Returns:
    #   tuple: hybridization enthalpy, enthropy and melting temperature.
    
    # Standard 1 M NaCl
    na_conc_0 = 1.

    # Calculate dH0(N-N)
    h = sum([tt[tt_mode][c][0] for c in couples])

    # Add initiation enthalpy
    if tt[tt_mode]['has_end']:
        h += tt[tt_mode]['end%s' % (seq[0],)][0]
        h += tt[tt_mode]['end%s' % (seq[-1],)][0]
    if tt[tt_mode]['has_init']:
        h += tt[tt_mode]['init'][0]

    # Correct enthalpy for symmetry
    if tt[tt_mode]['has_sym'] and seq == rc(seq, 'dna'):
        h += tt[tt_mode]['sym'][0]

    # Calculate dS0(N-N) in kcal / (K mol)
    s = sum([tt[tt_mode][c][1] for c in couples])

    # Add initiation enthropy
    if tt[tt_mode]['has_end']:
        s += tt[tt_mode]['end%s' % (seq[0],)][1]
        s += tt[tt_mode]['end%s' % (seq[-1],)][1]
    if tt[tt_mode]['has_init']:
        s += tt[tt_mode]['init'][1]

    # Correct enthalpy for symmetry
    if tt[tt_mode]['has_sym'] and seq == rc(seq, 'dna'):
        s += tt[tt_mode]['sym'][1]
    s /= 1e3

    # Calculate melting temperature in Celsius
    if 'sugimoto' in tt_mode:
        tm = h / (s + R * math.log(oligo_conc / 4))
    else:
        tm = h / (s + R * math.log(oligo_conc))

    # Output
    return((h, s, tm))

def duMelt_calc(name, seq, oligo_conc, na_conc, mg_conc,
	fa_conc, fa_mode, fa_mval_s, mvalue,
    tt, tt_mode, celsius, is_verbose,
    do_curve, curve_step, curve_range, curve_outpath):
    # Calculate melting temperature of provided oligo duplex sequence.
    # 
    # Args:
    #   seq (string): oligonucleotide sequence.
    #   oligo_conc (float): oligonucleotide concentration in M.
    #   na_conc (float): monovalent species concentration in M.
    #   mg_conc (float): divalent species concentration in M.
    #   fa_conc (float): formamide concentration in %v,v.
    #   fa_mode (string): formamide correction lavel.
    #   fa_mval_s (string): formamide m-value string.
    #   mvalue (lambda): formamide m-value function.
    #   tt (dict): thermodynamic table list.
    #   tt_mode (string): thermodynamic table label.
    #   celsius (bool): convert K to degC.
    #   is_verbose (bool): be verbose.
    #   do_curve (bool): generate melting curve.
    #   curve_step (float): melting curve temperature step.
    #   curve_range (float): melting curve temperature range.
    #   curve_outpath (string): melting curve output path.

    # Make string uppercase
    seq = seq.upper()

    # Check that all characters code for bases
    if not 0 == sum([1 for b in seq if b not in 'TCGAU']):
        print('The provided sequence contains non-nucleic acid characters.')
        return

    # Check that the correct nucleic acid type was provided
    if tt_mode in ['sugimotod', 'allawi']:
        if 'U' in seq:
            print('The option -t %s requires a DNA sequence.' % (tt_mode,))
            return
    elif tt_mode in ['sugimotor', 'freier']:
        if 'T' in seq:
            print('The option -t %s requires a RNA sequence.' % (tt_mode,))
            return

    # Make NN couples
    couples = [seq[i:(i+2)] for i in range(len(seq) - 1)]

    # Calculate GC content
    fgc = (seq.count('G') + seq.count('C')) / float(len(seq))

    # 1 M NaCl case
    # Based on SantaLucia, PNAS(95), 1998
    # -----------------------------------
    (h, s, Tm1) = duMelt_std_calc(seq, tt, tt_mode, couples, oligo_conc)

    # Adjust for FA
    # Based on Wright, Appl. env. microbiol.(80), 2014
    # Or on McConaughy, Biochemistry(8), 1969
    # ------------------------------------------------
    Tm2 = duMelt_fa_adj(Tm1, h, s, seq, oligo_conc,
    	fa_conc, fa_mode, mvalue, tt_mode)

    # Adjusted for [Na]
    # Based on Owczarzy et al, Biochemistry(43), 2004
    # -----------------------------------------------
    Tm3 = duMelt_na_adj(Tm2, na_conc, fgc)

    # Adjusted for Mg
    # Based on Owczarzy et al, Biochemistry(47), 2008
    # -----------------------------------------------
    if 0 < mg_conc:
        Tm4 = duMelt_mg_adj(Tm2, mg_conc, fgc)
    else:
        Tm4 = Tm3

    # Generate melting curves
    # -----------------------
    
    if do_curve:
        fout = open(curve_outpath, 'a+')
        tab = duMelt_curve(seq, oligo_conc, na_conc, mg_conc,
        	fa_conc, fa_mode, mvalue, fgc, h, s, Tm1,
        	curve_range, curve_step, tt_mode)
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
        print("%s\t%f\t%f\t%f\t%f\t%s" % (name, g, h, s, Tm4, seq))
    else:
        print("""
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
        ))

# CLASSES ======================================================================

# END ==========================================================================

################################################################################
