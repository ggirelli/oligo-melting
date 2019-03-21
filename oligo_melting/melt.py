
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for melting temperature calculation and correction for
              single-strand oligonucleotide secondary structures.
'''

import numpy as np
import oligo_melting as om
import os
import pandas as pd
import pkg_resources
import re
from tqdm import tqdm
import profile

R = 1.987 / 1000    # kcal / (K mol)
AB_DNA = ["ACGT", "TGCA"]
AB_RNA = ["ACGU", "UGCA"]
AB_NA = {"DNA":AB_DNA, "RNA":AB_RNA}

class Sequence(object):
    """docstring for Sequence"""
    def __init__(self, seq, t, name = None):
        super(Sequence, self).__init__()
        seq = seq.upper()
        self.text = seq
        self.len = len(seq)
        self.rc = self.rc(seq, t)
        self.fgc = (seq.count('G') + seq.count('C')) / self.len
        t = t.upper()
        assert t in AB_NA.keys(), "%s not in %s" % (t, list(AB_NA.keys()))
        self.natype = t
        self.ab = AB_NA[t]
        assert_msg = "sequence alphabet and nucleic acid type mismatch."
        assert_msg += "\n%s\t%s" % (set(seq), self.ab[0])
        assert all(x in self.ab[0] for x in set(seq)), assert_msg
        if type(name) == type(None):
            self.name = "%d-mer" % self.len
        else:
            self.name = name

    def dimers(self):
        '''Extract NN dimers from sequence.
        Args:
            seq (str): nucleic acid sequence.
        Returns:
            list: NN dimers.
        '''
        return (self.text[i:(i+2)] for i in range(self.len - 1))

    @staticmethod
    def rc(na, t):
        '''Calculate reverse complement.
        Args:
            na (string): nucleic acid sequence.
            t (string): nucleic acid type, either 'dna' or 'rna'.
        Return:
            string: reverse complement of na.
        '''
        t = t.upper()
        assert t in AB_NA.keys(), "%s not in %s" % (t, list(AB_NA.keys()))
        ab = AB_NA[t]

        rab = ab[1].strip().lower()
        ab = ab[0].strip().lower()
        na = na.lower()

        for c in na:
            assert_msg = 'provided string conflicts with selected alphabet.'
            assert_msg += " (%s in %s)" % (c, na)
            assert c in ab, assert_msg

        r = na[::-1]
        rc = []
        for c in r:
            rc.append(rab[ab.index(c)])
        rc = ''.join([str(c) for c in rc]).upper()

        return(rc)

    @staticmethod
    def dimerator(seq):
        '''Extract NN dimers from sequence.
        Args:
            seq (str): nucleic acid sequence.
        Returns:
            list: NN dimers.
        '''
        return (seq[i:(i+2)] for i in range(len(seq) - 1))

class NNEnergyTable(object):
    """docstring for NNEnergyTable"""
    def __init__(self, path, natypes):
        super(NNEnergyTable, self).__init__()

        natypes = natypes.upper()
        assert type(None) != type(re.match("[A-Z]+:[A-Z]+", natypes))
        assert all([x in NATYPES for x in natypes.split(":")])
        self.__natypes = natypes.split(":")

        assert os.path.isfile(path), "'%s' file not found." % path
        self.__table = pd.read_csv(path, "\t",
            index_col = 0, names = ["dH0", "dS0", "dG0"])

        self.has_init = "init" in self.__table.index
        if self.has_init:
            self.__init = self.__table.loc["init", :].to_dict()
            self.__table.drop("init", axis = 0)

        self.has_sym = "sym" in self.__table.index
        if self.has_sym:
            self.__sym = self.__table.loc["sym", :].to_dict()
            self.__table.drop("sym", axis = 0)

        self.has_end = any([x.startswith("end") for x in self.__table.index])
        if self.has_end:
            endRows = [x for x in self.__table.index if x.startswith("end")]
            self.__end = self.__table.loc[endRows, :]
            self.__table.drop(endRows, axis = 0)
            self.__end.index = [x[-1] for x in self.__end.index]
            self.__end = self.__end.to_dict()

        self.__ab = set("".join([x for x in self.__table.index if 2 == len(x)]))
        self.__table = self.__table.to_dict()

    @property
    def natypes(self):
        return self.__natypes

    @property
    def tt(self):
        return self.__table

    @property
    def dH0(self):
        return self.__table["dH0"]

    @property
    def dS0(self):
        return self.__table["dS0"]

    @property
    def dG0(self):
        return self.__table["dG0"]

    @property
    def end(self):
        if self.has_end:
            return self.__end.copy()

    @property
    def sym(self):
        if self.has_sym:
            return self.__sym.copy()

    @property
    def init(self):
        if self.has_init:
            return self.__init.copy()

    @property
    def ab(self):
        return self.__ab.copy()

NATYPES = ("DNA", "RNA")
NN_TABLES_PATH = {
    # Table from Freier et al, PNAS(83), 1986 - in 1 M NaCl [RNA]
    "RNA:RNA" : "nntables/freier.tsv",
    # Table from Sugimoto et al, Biochemistry(34), 1995 - in 1 M NaCl [DNA:RNA]
    # For calculation from the DNA sequence 5'-to-3'
    "DNA:RNA" : "nntables/sugimoto_2.tsv",
    # Table from Sugimoto et al, Biochemistry(34), 1995 - in 1 M NaCl [DNA:RNA]
    # For calculation from the RNA sequence 5'-to-3'
    "RNA:DNA" : "nntables/sugimoto_1.tsv",
    # Table from Allawi&Santalucia, Biochemistry(36), 1997 - in 1 M NaCl [DNA]
    "DNA:DNA" : "nntables/allawy.tsv"
}
PKGROOT = pkg_resources.resource_filename("oligo_melting", '')
NN_TABLES = [NNEnergyTable(os.path.join(PKGROOT, v), k)
    for k,v in NN_TABLES_PATH.items()]
NN_TABLES = dict([(":".join(x.natypes), x) for x in NN_TABLES])

class MeltingIonCorrector(object):
    """docstring for MeltingIonCorrector"""

    DEFAULT_MONOVALENT = 1
    DEFAULT_DIVALENT = 0
    __mono = DEFAULT_MONOVALENT
    __di = DEFAULT_DIVALENT

    def __init__(self):
        super(MeltingIonCorrector, self).__init__()

    @property
    def monovalent(self):
        return self.__mono
    @monovalent.setter
    def monovalent(self, conc):
        assert 0 <= conc
        self.__mono = conc

    @property
    def divalent(self):
        return self.__di
    @divalent.setter
    def divalent(self, conc):
        assert 0 <= conc
        self.__di = conc

    def correct_monovalent(self, tm, seq, conc0 = 1):
        '''Adjust melting temperature of a duplexx based on sodium concentration.
        Based on Owczarzy et al, Biochemistry(43), 2004
        Args:
          tm (float): melting temperature at [Na] = 1 M.
          seq (Sequence)
          conc0 (float): initial monovalent species concentration.
        Returns:
          float: adjusted melting temperature.
        '''
        if self.__mono == conc0:
            return tm

        a = 4.29e-5
        b = 3.95e-5
        c = 9.4e-6
        tm = 1./tm
        if 1 == self.__mono:
            tm -= (a*seq.fgc-b) * np.log(conc0/self.__mono)
            tm -= c*(np.log(conc0 / self.__mono)**2)
        elif 0 < self.__mono:
            tm += (a*seq.fgc-b) * np.log(self.__mono/conc0)
            tm += c*(np.log(self.__mono/conc0)**2)
        tm = 1./tm

        return tm

    def correct_divalent(self, tm, seq, correct = False):
        '''Adjust melting temperature a duplexx based on sodium concentration.
        Based on Owczarzy et al, Biochemistry(47), 2008
        Args:
          tm (float): melting temperature at [Mg] = 0 M and [Na] = 1 M.
          seq (Sequence)
        Returns:
          float: adjusted melting temperature.
        '''
        if self.__di == 0:
            return tm

        a = 3.92e-5
        b = -9.11e-6
        c = 6.26e-5 
        d = 1.42e-5
        e = -4.82e-4
        f = 5.25e-4
        g = 8.31e-5

        if correct:
            logConc = np.log(self.__mono)
            a *= 0.843-0.352*np.sqrt(self.__mono)*logConc
            d *= 1.279-4.03e-3*logConc-8.03e-3*logConc**2
            g *= 0.486-0.258*logConc+5.25e-3*logConc**3

        logConc = np.log(self.__di)
        tm = 1./tm + a + b*logConc + seq.fgc*(c + d*logConc)
        tm += (1./(2*(seq.len-1))) * (e + f*logConc + g*(logConc)**2)
        tm = 1./tm

        return tm

    def correct(self, tm, seq, mono_conc0 = 1):
        if 0 < self.__di:
            if 0 == self.__mono:
                return self.correct_divalent(tm, seq)
            else:
                ionRatio =  np.sqrt(self.__di)/self.__mono
                if ionRatio < 0.22:
                    return self.correct_monovalent(tm, seq, mono_conc0)
                elif ionRatio < 6:
                    return self.correct_divalent(tm, seq, True)
                else:
                    return self.correct_divalent(tm, seq)
        elif 0 < self.__mono:
            return self.correct_monovalent(tm, seq, mono_conc0)
        else:
            return tm

class MeltingDenaturantCorrector(object):
    """docstring for MeltingDenaturantCorrector"""

    MODES = ("MCCONAUGHY", "WRIGHT")
    DEFAULT_MODE = MODES[0]
    DEFAULT_CONC = 0
    DEFAULT_M1 = 0.1734
    DEFAULT_M2 = 0
    __denaturant = DEFAULT_CONC
    __m1 = DEFAULT_M1
    __m2 = DEFAULT_M2
    __mode = DEFAULT_MODE

    def __init__(self):
        super(MeltingDenaturantCorrector, self).__init__()

    @property
    def conc(self):
        return self.__denaturant
    @conc.setter
    def conc(self, conc):
        assert 0 <= conc and 100 >= conc
        self.__denaturant = conc

    @property
    def m1(self):
        return self.m1
    @m1.setter
    def m1(self, m):
        self.__m1 = m

    @property
    def m2(self):
        return self.__denaturant
    @m2.setter
    def m2(self, m):
        self.__m2 = m

    @property
    def mode(self):
        return self.__mode
    @mode.setter
    def mode(self, mode):
        assert mode in self.MODES
        self.__mode = mode
    

    def mcconaughy_correction(self, tm, conc0 = 0, **kwargs):
        '''Adjust melting temperature of a duplex based on formamide
        concentration. Based on McConaughy, Biochemistry(8), 1969
        Args:
          tm (float): melting temperature.
          conc0 (float): initial formamide concentration in %v,v.
        Returns:
          float: corrected melting temperature.
        '''
        deltaDenaturant = self.__denaturant - conc0
        if 0 == deltaDenaturant:
            return tm
        else:
            tm -= 0.72 * deltaDenaturant
        return tm

    def wright_correction(self, tm, h, s, seq, oligo, conc0 = 0, **kwargs):
        '''Adjust melting temperature of a duplex based on formamide
        concentration. Based on Wright, Appl. env. microbiol.(80), 2014
        Args:
          tm (float): melting temperature.
          h (float): standard enthalpy, in kcal / mol.
          s (float): standard enthropy, in kcal / (K mol).
          seq (Sequence)
          oligo (float): oligonucleotide concentration in M.
          conc0 (float): initial formamide concentration in %v,v.
        Returns:
          float: corrected melting temperature.
        '''
        deltaDenaturant = self.__denaturant - conc0
        if 0 == deltaDenaturant:
            return tm

        tm = (h + self.mvalue() * deltaDenaturant) / (R * np.log(oligo) + s)
        return tm

    def mvalue(self):
        return self.__m1 if 0 == self.__m2 else self.__m1 * seq.len + self.__m2

    def correct(self, tm, h, s, seq, oligo, conc0 = 0):
        return getattr(self, "%s_correction" % self.__mode.lower())(
            seq = seq, h = h, s = s, tm = tm, oligo = oligo, conc0 = conc0
        )

class Melter(object):
    """docstring for Melter"""

    DEGREES_TYPE = ["CELSIUS", "KELVIN"]
    DEFAULT_DEGREES_TYPE = DEGREES_TYPE[1]
    DEFAULT_NN = 'DNA:DNA'
    DEFAULT_OLIGO = .25e-6
    DEFAULT_RANGE = 10
    DEFAULT_STEP = .1

    __oligo = DEFAULT_OLIGO
    ions = MeltingIonCorrector()
    denaturant = MeltingDenaturantCorrector()
    __nnet = NN_TABLES[DEFAULT_NN]
    __degrees = DEFAULT_DEGREES_TYPE
    __curve = [DEFAULT_RANGE, DEFAULT_STEP]

    def __init__(self):
        super(Melter, self).__init__()

    @property
    def nnet(self):
        return self.__nnet
    @nnet.setter
    def nnet(self, tt):
        assert NNEnergyTable == type(tt)
        self.__nnet = tt

    @property
    def oligo(self):
        return self.__oligo
    @oligo.setter
    def oligo(self, conc):
        assert 0 <= conc
        if self.__nnet.natypes[0] != self.__nnet.natypes[1]:
            conc /= 4
        self.__oligo = conc

    @property
    def degrees(self):
        return self.__degrees
    @degrees.setter
    def degrees(self, d):
        assert d in self.DEGREES_TYPE
        self.__degrees = d
    
    @property
    def curve_range(self):
        return self.__curve[0]
    @curve_range.setter
    def curve_range(self, x):
        self.__curve[0] = x

    @property
    def curve_step(self):
        return self.__curve[1]
    @curve_step.setter
    def curve_step(self, x):
        self.__curve[1] = x

    def load_nn_table(self, nntype):
        assert nntype in NN_TABLES.keys()
        self.__nnet = NN_TABLES[nntype]

    def calculate(self, seq):
        '''Calculate melting temperature of provided sequence.
        Args:
            seq (Sequence/string)
        Returns:
            tuple: sequence name, hybridization enthalpy, enthropy, melting
                   temperature, and sequence
        '''
        if not type(seq) == Sequence:
            seq = Sequence(seq, self.nnet.natypes[0])

        # 1 M NaCl case; SantaLucia, PNAS(95), 1998
        name, g, h, s, tmStd, text = self.__calculate_standard(seq, True)

        # Adjust for FA; Wright, Appl. env. microbiol.(80), 2014
        # Or McConaughy, Biochemistry(8), 1969
        tm = self.denaturant.correct(tmStd, h, s, seq, self.oligo)

        # Adjust for [Na]; Owczarzy et al, Biochemistry(43), 2004
        # Adjust for Mg; Owczarzy et al, Biochemistry(47), 2008
        tm = self.ions.correct(tm, seq)

        g = h - (37 + 273.15) * s
        if self.degrees == self.DEGREES_TYPE[0]:
            tm -= 273.15
        return((name, g, h, s, tm, text))

    def __calculate_standard(self, seq, forceKelvin = False):
        '''Calculate melting temperature of a duplex at standard 1 M NaCl
        (monovalent ions conc). Based on SantaLucia, PNAS(95), 1998
        Args:
            seq (Sequence/string)
        Returns:
            tuple: sequence name, hybridization enthalpy, entropy, melting
                   temperature, and sequence
        '''

        h = sum([self.__nnet.dH0[c] for c in seq.dimers()])
        s = sum([self.__nnet.dS0[c] for c in seq.dimers()])
        if self.__nnet.has_end:
            h += self.__nnet.end["dH0"][seq.text[0]]
            h += self.__nnet.end["dH0"][seq.text[-1]]
            s += self.__nnet.end["dS0"][seq.text[0]]
            s += self.__nnet.end["dS0"][seq.text[-1]]
        if self.__nnet.has_init:
            h += self.__nnet.init["dH0"]
            s += self.__nnet.init["dS0"]
        if self.__nnet.has_sym and seq == seq.rc:
            h += self.__nnet.sym["dH0"]
            s += self.__nnet.sym["dS0"]
        s /= 1e3

        tm = h / (s + R * np.log(self.__oligo))

        g = h - (37 + 273.15) * s
        if self.degrees == self.DEGREES_TYPE[0] and not forceKelvin:
            tm -= 273.15
        return((seq.name, g, h, s, tm, seq.text))

    def melting_curve(self, seq):
        '''Generate melting curve.
        Args:
            seq (Sequence)
            g (float): hybridization free energy
            h (float): hybridization enthalpy
            s (float): hybridization entropy
            tm (float): melting temperature
        Returns:
            list of tuples: [(temperature, dissociated fraction), ...]
        '''

        if not type(seq) == Sequence:
            seq = Sequence(seq, self.nnet.natypes[0])
        
        name, g, h, s, tm, text = self.__calculate_standard(seq, True)

        if self.denaturant.mode == self.denaturant.MODES[1]:
            tm = self.denaturant.correct(tm, h, s, seq, self.oligo)
            def compute_curve_step(t, h, s, seq):
                m = self.denaturant.mvalue()
                k = self.__dissoc_fraction(t, h, s, m * self.denaturant.conc)
                return (t, k)
        else:
            def compute_curve_step(t, h, s, seq):
                k = self.__dissoc_fraction(t, h, s, 0)
                t = self.denaturant.correct(t, h, s, seq, self.oligo)
                return (t, k)

        if self.degrees == self.DEGREES_TYPE[0]:
            parse_tm = lambda x: x - 273.15
        else:
            parse_tm = lambda x: x

        curve = []
        tstart = tm - self.curve_range / 2
        for ti in range(int(self.curve_range/self.curve_step)):
            t, f = compute_curve_step(tstart + self.curve_step * ti, h, s, seq)
            t = self.ions.correct(t, seq)
            curve.append((parse_tm(t), f))
        
        curve = np.array(curve)

        return curve

    def __dissoc_fraction(self, t, h, s, gplus):
        '''Calculate duplex dissociated fraction at given temperature.
        Args:
            t (float): temperature
            h (float): enthalpy
            s (float): entropy
            g (float): free energy difference due to denaturants
        Returns:
            float: dissociated fraction
        '''
        dg = h - t * s + gplus
        factor = np.exp(-dg / (R * t)) * self.oligo
        return 1 / (1 + factor)
