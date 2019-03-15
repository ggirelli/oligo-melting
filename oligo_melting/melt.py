
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
import re

R = 1.987 / 1000    # kcal / (K mol)
CELSIUS = "Celsius"
KELVIN = "Kelvin"
AB_DNA = ["ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN"]
AB_RNA = ["ACGURYKMSWBDHVN", "UGCAYRMKSWVHDBN"]
AB_NA = {"DNA":AB_DNA, "RNA":AB_RNA}

class Sequence(object):
    """docstring for Sequence"""
    def __init__(self, seq, t):
        super(Sequence, self).__init__()
        seq = seq.upper()
        self.text = seq
        self.len = len(seq)
        self.rc = self.rc(seq, t)
        self.fgc = (seq.count('G') + seq.count('C')) / float(len(seq))
        t = t.upper()
        assert t in AB_NA.keys(), "%s not in %s" % (t, list(AB_NA.keys()))
        self.natype = t
        self.ab = AB_NA[t]
        assert_msg = "sequence alphabet and nucleic acid type mismatch."
        assert all(x in self.ab[0] for x in set(seq)), assert_msg

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
            assert c in ab, 'provided string conflicts with selected alphabet.'

        r = na[::-1]
        rc = []
        for c in r:
            rc.append(rab[ab.index(c)])
        rc = ''.join([str(c) for c in rc]).upper()

        return(rc)

    @staticmethod
    def get_dimers(seq):
        '''Extract NN dimers from sequence.
        Args:
            seq (str): nucleic acid sequence.
        Returns:
            list: NN dimers.
        '''
        return([seq[i:(i+2)] for i in range(len(seq) - 1)])

class NNEnergyTable(object):
    """docstring for NNEnergyTable"""
    def __init__(self, path, natypes):
        super(NNEnergyTable, self).__init__()

        natypes = natypes.upper()
        assert type(None) != type(re.match("[A-Z]+:[A-Z]+", natypes))
        assert all([x in NATYPES for x in natypes.split(":")])
        self.__natypes = natypes.split(":")

        assert os.path.isfile(path)
        self.__table = pd.read_csv(path, "\t",
            index_col = 0, names = ["dH0", "dS0", "dG0"])

        self.has_init = "init" in self.__table.index
        if self.has_init:
            self.__init = self.__table.loc["init", :]
            self.__table.drop("init", axis = 0)

        self.has_sym = "sym" in self.__table.index
        if self.has_sym:
            self.__sym = self.__table.loc["sym", :]
            self.__table.drop("sym", axis = 0)

        self.has_end = any([x.startswith("end") for x in self.__table.index])
        if self.has_end:
            endRows = [x for x in self.__table.index if x.startswith("end")]
            self.__end = self.__table.loc[endRows, :]
            self.__table.drop(endRows, axis = 0)
            self.__end.index = [x[-1] for x in self.__end.index]

        self.__ab = set("".join([x for x in self.__table.index if 2 == len(x)]))

    @property
    def natypes(self):
        return self.__natypes

    @property
    def tt(self):
        return self.__table.copy()

    @property
    def dH0(self):
        return self.__table.copy().loc[:, "dH0"]

    @property
    def dS0(self):
        return self.__table.copy().loc[:, "dS0"]

    @property
    def dG0(self):
        return self.__table.copy().loc[:, "dG0"]

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
NN_TABLES = [NNEnergyTable(v, k) for k,v in NN_TABLES_PATH.items()]
NN_TABLES = dict([(":".join(x.natypes), x) for x in NN_TABLES])

class MeltingIonCorrector(object):
    """docstring for MeltingIonCorrector"""

    __mono = .3
    __di = 0

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

    __denaturant = 25

    def __init__(self):
        super(MeltingDenaturantCorrector, self).__init__()

    @property
    def conc(self):
        return self.__denaturant
    @conc.setter
    def conc(self, conc):
        assert 0 <= conc and 100 >= conc
        self.__denaturant = conc

    def wright_correction(self, tm, conc0 = 0):
        '''Adjust melting temperature of a duplex based on formamide concentration.
        Based on Wright, Appl. env. microbiol.(80), 2014
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

    def mcconaughy_correction(self, h, s, m, oligo, conc0 = 0):
        '''Adjust melting temperature of a duplex based on formamide concentration.
        Based on McConaughy, Biochemistry(8), 1969
        Args:
          h (float): standard enthalpy, in kcal / mol.
          s (float): standard enthropy, in kcal / (K mol).
          oligo (float): oligonucleotide concentration in M.
          m (float): formamide m-value function.
          conc0 (float): initial formamide concentration in %v,v.
        Returns:
          float: corrected melting temperature.
        '''
        deltaDenaturant = self.__denaturant - conc0
        if 0 == deltaDenaturant:
            return tm
        else:
            if self.__nnet.natypes[0] != self.__nnet.natypes[1]:
                oligo /= 4
            tm = (h + m * deltaDenaturant) / (R * np.log(oligo) + s)
        return tm

class Melter(object):
    """docstring for Melter"""

    __oligo = .25e-6
    ions = MeltingIonCorrector()
    denaturant = MeltingDenaturantCorrector()
    __nnet = NN_TABLES['DNA:DNA']
    __degrees = "Celsius" # Either 'Celsius' or 'Kelvin'
    __verbose = False

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
        self.__oligo = conc
    
    def calculate(self, seq):
        seq = Sequence(seq, self.nnet.natypes[0])

        # 1 M NaCl case; SantaLucia, PNAS(95), 1998
        h0, s0, tmStd = self.__calculate_standard(seq, True)
        print(("tmStd", tmStd -273.15))

        # Adjust for FA; Wright, Appl. env. microbiol.(80), 2014
        tm = self.denaturant.wright_correction(tmStd)
        print(("tm", tm -273.15))
        # Or McConaughy, Biochemistry(8), 1969 !!!!!!!!!!!!!!!!! <- to implement

        # Adjust for [Na]; Owczarzy et al, Biochemistry(43), 2004
        # Adjust for Mg; Owczarzy et al, Biochemistry(47), 2008
        tm = self.ions.correct(tm, seq)
        print(("tm", tm -273.15))

        if self.__degrees == CELSIUS:
            tm -= 273.15
        return(tm)

    def __calculate_standard(self, seq, forceKelvin = False):
        '''Calculate melting temperature of a duplex at standard 1 M NaCl
        (monovalent ions conc). Based on SantaLucia, PNAS(95), 1998
        Args:
          seq (string): oligonucleotide sequence.
        Returns:
          tuple: hybridization enthalpy, enthropy and melting temperature.'''
        dimers = seq.get_dimers(seq.text)
        monovalent = 1 # 1 M

        h = sum([self.__nnet.dH0[c] for c in dimers])
        s = sum([self.__nnet.dS0[c] for c in dimers])
        if self.__nnet.has_end:
            h += self.__nnet.end.loc[[seq.text[0], seq.text[-1]], "dH0"].sum()
            s += self.__nnet.end.loc[[seq.text[0], seq.text[-1]], "dS0"].sum()
        if self.__nnet.has_init:
            h += self.__nnet.init["dH0"]
            s += self.__nnet.init["dS0"]
        if self.__nnet.has_sym and seq == seq.rc:
            h += self.__nnet.sym["dH0"]
            s += self.__nnet.sym["dS0"]
        s /= 1e3

        if self.__nnet.natypes[0] != self.__nnet.natypes[1]:
            self.__oligo /= 4
        tm = h / (s + R * np.log(self.__oligo))

        if self.__degrees == CELSIUS and not forceKelvin:
            tm -= 273.15

        return((h, s, tm))

m = Melter()
m.ions.monovalent = 0.3
m.ions.divalent = 0
m.denaturant.conc = 0
#m.nnet = NN_TABLES["DNA:RNA"]
print(m.calculate("ATGTCAATGTCAATGTCAATGTCAATGTCAATGTCAATGTCAATGTCA"))
