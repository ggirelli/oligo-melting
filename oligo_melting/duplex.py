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

# FUNCTIONS ====================================================================

def melt_fa_adj(tm, h, s, seq, oligo_conc,
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
    if OligoMelt.FA_MODE_MCCONA1969 == fa_mode:
        tm -= 0.72 * dfac

    # Apply to free energy
    if OligoMelt.FA_MODE_WRIGHT2014 == fa_mode:
        m = mvalue(len(seq))
        if tt_mode in OligoMelt.NN_HYB_LABELS:
            tm = (h + m * dfac) / (R * math.log(oligo_conc / 4) + s)
        else:
            tm = (h + m * dfac) / (R * math.log(oligo_conc) + s)

    return(tm)

def melt_na_adj(tm, na_conc, fgc, na_conc_0 = None):
    # Adjust melting temperature of a duplexx based on sodium concentration.
    # Based on Owczarzy et al, Biochemistry(43), 2004
    # 
    # Args:
    #   tm (float): melting temperature at [Na] = 1 M.
    #   na_conc (float): monovalent species concentration in M.
    #   fgc (float): GC content fraction.
    #   na_conc_0 (float): initial monovalent species concentration.
    #   
    # Returns:
    #   float: adjusted melting temperature.
    
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

def melt_mg_adj(tm, mg_conc, fgc):
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

def melt_ion_adj(tm, na_conc, mg_conc, fgc, na_conc_0 = None):
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
    
    if type(None) == type(na_conc_0):
        na_conc_0 = 1.

    if 0 != mg_conc:
        return(melt_mg_adj(tm, mg_conc, fgc))
    elif 0 != na_conc:
        return(melt_na_adj(tm, na_conc, fgc, na_conc_0))
    else:
        return(tm)

def melt_curve(seq, oligo_conc, na_conc, mg_conc, fa_conc, fa_mode, mvalue,
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
    if OligoMelt.FA_MODE_WRIGHT2014 == fa_mode:
        tm = melt_fa_adj(tm, h, s, seq, oligo_conc,
        	fa_conc, fa_mode, mvalue, tt_mode)

    # Explore the temperature range
    t = tm - trange / 2.
    while t <= tm + trange / 2.:
        if OligoMelt.FA_MODE_MCCONA1969 == fa_mode:
            # Calculate dissociated fraction
            k = dissoc_fraction(h, t, s, oligo_conc, tt_mode, 0)

            # Adjust output temperature
            t_out = melt_fa_adj(t, h, s, seq, oligo_conc,
            	fa_conc, fa_mode, mvalue, tt_mode)
            t_out = melt_ion_adj(t_out, na_conc, mg_conc, fgc)

        if OligoMelt.FA_MODE_WRIGHT2014 == fa_mode:
            # Calculate current FA m-value
            m = mvalue(len(seq))

            # Calculate dissociated fraction
            k = dissoc_fraction(h, t, s, oligo_conc, tt_mode, m * fa_conc)

            # Adjust output temperature
            t_out = melt_ion_adj(t, na_conc, mg_conc, fgc)

        # Append melting data
        data.append((t_out, k))
        
        # Go to next temperature
        t += tstep

    # Return melting table
    return(data)

def dissoc_fraction(h, t, s, oligo_conc, tt_mode, gplus):
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
    if tt_mode in OligoMelt.NN_HYB_LABELS:
        factor = math.exp(-dg / (R * t)) * (oligo_conc / 4)
    else:
        factor = math.exp(-dg / (R * t)) * oligo_conc

    # Calculate fraction
    k = 1 - factor / (1 + factor)

    # Output
    return(k)

def melt_std_calc(seq, tt, tt_mode, couples, oligo_conc):
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
    if tt[tt_mode]['has_sym'] and seq == OligoMelt.rc(seq, 'dna'):
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
    if tt[tt_mode]['has_sym'] and seq == OligoMelt.rc(seq, 'dna'):
        s += tt[tt_mode]['sym'][1]
    s /= 1e3

    # Calculate melting temperature in Celsius
    if tt_mode in OligoMelt.NN_HYB_LABELS:
        tm = h / (s + R * math.log(oligo_conc / 4))
    else:
        tm = h / (s + R * math.log(oligo_conc))

    # Output
    return((h, s, tm))

def melt_calc(name, seq, oligo_conc, na_conc, mg_conc,
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
    if tt_mode in OligoMelt.NN_DNA_TEMPLATE_LABELS:
        if 'U' in seq:
            print('The option -t %s requires a DNA sequence.' % (tt_mode,))
            return
    elif tt_mode in OligoMelt.NN_RNA_TEMPLATE_LABELS:
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
    (h, s, Tm1) = melt_std_calc(seq, tt, tt_mode, couples, oligo_conc)

    # Adjust for FA
    # Based on Wright, Appl. env. microbiol.(80), 2014
    # Or on McConaughy, Biochemistry(8), 1969
    # ------------------------------------------------
    Tm2 = melt_fa_adj(Tm1, h, s, seq, oligo_conc,
    	fa_conc, fa_mode, mvalue, tt_mode)

    # Adjusted for [Na]
    # Based on Owczarzy et al, Biochemistry(43), 2004
    # -----------------------------------------------
    Tm3 = melt_na_adj(Tm2, na_conc, fgc)

    # Adjusted for Mg
    # Based on Owczarzy et al, Biochemistry(47), 2008
    # -----------------------------------------------
    if 0 < mg_conc:
        Tm4 = melt_mg_adj(Tm2, mg_conc, fgc)
    else:
        Tm4 = Tm3

    # Generate melting curves
    # -----------------------
    
    if do_curve:
        fout = open(curve_outpath, 'a+')
        tab = melt_curve(seq, oligo_conc, na_conc, mg_conc,
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

# END ==========================================================================

################################################################################
