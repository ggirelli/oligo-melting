# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Description: keeps constant values for the oligo-melting package.
# 
# ------------------------------------------------------------------------------

# DEPENDENCIES =================================================================

import pkg_resources

# FUNCTIONS ====================================================================

# Taken from:
# http://code.activestate.com/recipes/65207-constants-in-python/?in=user-97991
class _const:
	class ConstError(TypeError): pass
	def __setattr__(self, name, value):
		if name in self.__dict__.keys():
			raise(self.ConstError, "Can't rebind const(%s)" % name)
		self.__dict__[name] = value

# CONSTANTS ====================================================================

# Alphabets --------------------------------------------------------------------

_const.AB_DNA = ["ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN"]
_const.AB_RNA = ["ACGURYKMSWBDHVN", "UGCAYRMKSWVHDBN"]

# Constants --------------------------------------------------------------------

# Gas constant
_const.R = 1.987 / 1000    # kcal / (K mol)

# Thermodynamic tables ---------------------------------------------------------

# Table from Freier et al, PNAS(83), 1986 - in 1 M NaCl [RNA]
# dH0: kcal / mol
# dS0: eu = cal / (K mol)
# dG0: kcal / mol
_const.NN_TABLE_RNA_RNA = {
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
_const.NN_TABLE_DNA_RNA = {
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
_const.NN_TABLE_RNA_DNA = {
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
_const.NN_TABLE_DNA_DNA = {
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
_const.NN_TABLES = {
    "RNA:RNA" : _const.NN_TABLE_RNA_RNA, "DNA:DNA" : _const.NN_TABLE_DNA_DNA,
    "DNA:RNA" : _const.NN_TABLE_DNA_RNA, "RNA:DNA" : _const.NN_TABLE_RNA_DNA
}
_const.NN_LABELS = list(_const.NN_TABLES.keys())
_const.NN_HYB_LABELS = ["DNA:RNA", "RNA:DNA"]
_const.NN_DNA_TEMPLATE_LABELS = ["DNA:RNA", "DNA:DNA"]
_const.NN_RNA_TEMPLATE_LABELS = ["RNA:RNA", "RNA:DNA"]

# Formamide correction ---------------------------------------------------------

_const.FA_MODE_LABELS = ["wright", "mcconaughy"]
_const.FA_MODE_WRIGHT2014 = _const.FA_MODE_LABELS[0]
_const.FA_MODE_MCCONA1969 = _const.FA_MODE_LABELS[1]

# RUN ==========================================================================

# Save constants
import sys
sys.modules[__name__] = _const()

# END ==========================================================================

################################################################################
