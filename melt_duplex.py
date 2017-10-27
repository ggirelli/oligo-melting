#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.4.4
# Date: 20170711
# Project: oligo characterization
# Description:  calculate melting temperature of a provide NA duplex
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

import argparse
import os
import re
import sys

from lib.meltlib import *

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(
    description = '''
Calculate melting temperature of a DNA duplex at provided [oligo],
[Na+], [Mg2+]. Either provide an oligo sequence or a file with one oligo
per line (and use -F option). References:
 [1] Freier et al, PNAS(83), 1986;
 [2] Sugimoto et al, Biochemistry(34), 1995.
 [3] Allawi & Santalucia, Biochemistry(36), 1997;
 [4] SantaLucia, PNAS(95), 1998;
 [5] Owczarzy et al, Biochemistry(43), 2004;
 [6] Owczarzy et al, Biochemistry(47), 2008;
 [7] McConaughy et al, Biochemistry(8), 1969;
 [8] Wright et al, Appl. env. microbiol.(80), 2014.
''')

# Add mandatory arguments
parser.add_argument('seq', type = str, nargs = 1, help = '''
    DNA duplex sequence (one strand only) or path to a FASTA file (use with -F).
    ''')

# Add arguments with default value
parser.add_argument('-t', '--type', type = str, nargs = 1,
    help = '''Duplex type. Possible values: DNA:DNA (based on ref.3, default),
    RNA:RNA (based on ref.1), DNA:RNA (based on ref.2., given DNA sequence)
    or RNA:DNA (based on ref.2, given RNA sequence). The first nucleic acid type
    indicates the provided sequence.''',
    choices = ['DNA:DNA', 'RNA:RNA', 'RNA:DNA', 'DNA:RNA'],
    default = ['DNA:DNA'])
parser.add_argument('-o', '--oconc', metavar = "oligo_conc",
    type = float, nargs = 1,
    help = '''Oligonucleotide concentration [M].
    Default: 0.25e-6 M''',
    default = [0.25e-6])
parser.add_argument('-n', '--naconc', metavar = "na_conc",
    type = float, nargs = 1,
    help = '''Na+ concentration [M].
    Default: 50e-3 M''',
    default = [50e-3])
parser.add_argument('-m', '--mgconc', metavar = "mg_conc",
    type = float, nargs = 1,
    help = '''Mg2+ concentration [M]. Note: Mg2+ correction overwrites Na+
    correction. Default: 0 M''',
    default = [0])
parser.add_argument('-f', '--faconc', type = float, nargs = 1,
    metavar = 'fa_conc', help = '''Formamide concentration in %%(v,v).''',
    default = [0])
parser.add_argument('--fa-mode', type = str, nargs = 1,
    metavar = 'fa_mode', help = '''Mode of formamide correction. "mcconaughy"
    for classical -0.72%%FA correction from ref. 7, "wright" for single reaction
    model correction from ref.8 (default).''',
    choices = ['mcconaughy', 'wright'], default = ['wright'])
parser.add_argument('--fa-mvalue', type = str, nargs = 1, metavar = 'm',
    help = '''Specify the formamide m-value to be used with the wright
    correction model. Use either a single value "x" or two values "xL+y" where
    L is the probe length. Default: 0.1734''', default = ["0.1734"])
parser.add_argument('--t-curve', type = float, nargs = 2,
    metavar = ('range', 'step'), help = '''Temperature range and step for
    melting curve generation. Use --make-curve option to generate the curve.
    Default: 20 degC range and 0.5 degC step.''', default = [20, 0.5])
parser.add_argument('--out-curve', type = str, nargs = 1, metavar = "outname",
    help = '''Path to output table containing tabulated curves.''',
    default = [False])

# Add flags
parser.add_argument('-C', '--celsius',
    dest = 'celsius', action = 'store_const',
    const = True, default = False,
    help = 'Output temperature in Celsius degrees. Default: Kelvin')
parser.add_argument('-F', '--usefile',
    dest = 'usefile', action = 'store_const',
    const = True, default = False,
    help = 'Use when a file path is provided instead of a single sequence.')
parser.add_argument('-v', '--verbose',
    dest = 'verbose', action = 'store_const',
    const = True, default = False,
    help = 'Verbose output.')

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables ------------------------------------------------

# Input oligo
seq = args.seq[0]

# Concentrations
oligo_conc = args.oconc[0]
na_conc = args.naconc[0]
mg_conc = args.mgconc[0]

# Thermodynamic table
tt_mode = args.type[0]

# Temperature units
celsius = args.celsius

# File as input
use_file = args.usefile

# Verbose mode
is_verbose = args.verbose

# Melting curve
curve_range = args.t_curve[0]
curve_step = args.t_curve[1]
curve_outpath = args.out_curve[0]
do_curve = False != curve_outpath

# Formamide
fa_conc = args.faconc[0]
fa_mode = args.fa_mode[0]
fa_mval_s = args.fa_mvalue[0]

# Additional checks ------------------------------------------------------------

# Check proper curve step/range pair
if curve_step > curve_range / 2:
    sys.exit("!!!ERROR! Curve step must be smaller than curve range.")
curve_range -= curve_range % curve_step

# Remove output curve file if it exists
if do_curve and os.path.isfile(curve_outpath):
    os.remove(curve_outpath)

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
    sys.exit("!!!ERROR! Unexpected formamide m-value format Check help page.")

# Fix m-value label
if not fa_mode in ["wright"]:
    fa_mval_s = "-"

# Thermodynamic tables
tt = {
    'RNA:RNA'   : freier,
    'DNA:DNA'   : allawi,
    'RNA:DNA'   : sugimotor,
    'DNA:RNA'   : sugimotod
}

# FUNCTIONS ====================================================================

# RUN ==========================================================================

# Build argument dictionary
data = {
    'oligo_conc' : oligo_conc,
    'na_conc' : na_conc,
    'mg_conc' : mg_conc,
    'fa_conc' : fa_conc,
    'fa_mode' : fa_mode,
    'fa_mval_s' : fa_mval_s,
    'mvalue' : mvalue,
    'tt' : tt,
    'tt_mode' : tt_mode,
    'celsius' : celsius,
    'is_verbose' : is_verbose,
    'do_curve' : do_curve,
    'curve_step' : curve_step,
    'curve_range' : curve_range,
    'curve_outpath' : curve_outpath
}

# CALCULATE --------------------------------------------------------------------

if not use_file:
    # Single sequence case
    data['name'] = 'seq'
    data['seq'] = seq
    duMelt_calc(**data)
else:
    # Empty sequence dictionary
    fastad = {}

    # Input file case
    curr_head = ""
    with open(seq) as fin:
        for row in fin:
            if ">" == row[0]:
                curr_head = row[1:].strip()
            else:
                if curr_head in fastad.keys():
                    fastad[curr_head] = fastad[curr_head] + row.strip()
                else:
                    fastad[curr_head] = row.strip()

    # Calculate for each fasta item
    for (name, seq) in fastad.items():
        data['name'] = name
        data['seq'] = seq
        duMelt_calc(**data)

# END ==========================================================================

################################################################################
