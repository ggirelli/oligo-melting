"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for melting temperature calculation and correction for
              single-strand oligonucleotide secondary structures.
"""

import numpy as np  # type: ignore
from oligo_melting import const
from oligo_melting.entab import NN_TABLES, NNEnergyTable
from oligo_melting.sequence import Sequence  # type: ignore


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

    def correct_monovalent(self, tm, seq, conc0=1):
        """Adjust melting temperature of a duplexx based on sodium concentration.
        Based on Owczarzy et al, Biochemistry(43), 2004
        Args:
          tm (float): melting temperature at [Na] = 1 M.
          seq (Sequence)
          conc0 (float): initial monovalent species concentration.
        Returns:
          float: adjusted melting temperature.
        """
        if self.__mono == conc0:
            return tm

        a = 4.29e-5
        b = 3.95e-5
        c = 9.4e-6
        tm = 1.0 / tm
        if 1 == self.__mono:
            tm -= (a * seq.fgc - b) * np.log(conc0 / self.__mono)
            tm -= c * (np.log(conc0 / self.__mono) ** 2)
        elif 0 < self.__mono:
            tm += (a * seq.fgc - b) * np.log(self.__mono / conc0)
            tm += c * (np.log(self.__mono / conc0) ** 2)
        tm = 1.0 / tm

        return tm

    def correct_divalent(self, tm, seq, correct=False):
        """Adjust melting temperature a duplexx based on sodium concentration.
        Based on Owczarzy et al, Biochemistry(47), 2008
        Args:
          tm (float): melting temperature at [Mg] = 0 M and [Na] = 1 M.
          seq (Sequence)
        Returns:
          float: adjusted melting temperature.
        """
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
            a *= 0.843 - 0.352 * np.sqrt(self.__mono) * logConc
            d *= 1.279 - 4.03e-3 * logConc - 8.03e-3 * logConc ** 2
            g *= 0.486 - 0.258 * logConc + 5.25e-3 * logConc ** 3

        logConc = np.log(self.__di)
        tm = 1.0 / tm + a + b * logConc + seq.fgc * (c + d * logConc)
        tm += (1.0 / (2 * (seq.len - 1))) * (e + f * logConc + g * (logConc) ** 2)
        tm = 1.0 / tm

        return tm

    def __correct_multivalent(self, tm, seq, mono_conc0=1):
        ionRatio = np.sqrt(self.__di) / self.__mono
        if ionRatio < 0.22:
            return self.correct_monovalent(tm, seq, mono_conc0)
        elif ionRatio < 6:
            return self.correct_divalent(tm, seq, True)
        else:
            return self.correct_divalent(tm, seq)

    def correct(self, tm, seq, mono_conc0=1):
        if 0 < self.__di:
            if 0 == self.__mono:
                return self.correct_divalent(tm, seq)
            else:
                return self.__correct_multivalent(tm, seq, mono_conc0)
        elif 0 < self.__mono:
            return self.correct_monovalent(tm, seq, mono_conc0)
        else:
            return tm


class MeltingDenaturantCorrector(object):
    """docstring for MeltingDenaturantCorrector"""

    DEFAULT_MODE = const.DENATURANT_MODES.MCCONAUGHY
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
        assert mode in const.DENATURANT_MODES
        self.__mode = mode

    def mcconaughy_correction(self, tm, conc0=0, **kwargs):
        """Adjust melting temperature of a duplex based on formamide
        concentration. Based on McConaughy, Biochemistry(8), 1969
        Args:
          tm (float): melting temperature.
          conc0 (float): initial formamide concentration in %v,v.
        Returns:
          float: corrected melting temperature.
        """
        deltaDenaturant = self.__denaturant - conc0
        if 0 == deltaDenaturant:
            return tm
        else:
            tm -= 0.72 * deltaDenaturant
        return tm

    def wright_correction(self, tm, h, s, seq, oligo, conc0=0, **kwargs):
        """Adjust melting temperature of a duplex based on formamide
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
        """
        deltaDenaturant = self.__denaturant - conc0
        if 0 == deltaDenaturant:
            return tm

        tm = (h + self.mvalue(seq) * deltaDenaturant) / (const.R * np.log(oligo) + s)
        return tm

    def mvalue(self, seq):
        return self.__m1 if 0 == self.__m2 else self.__m1 * seq.len + self.__m2

    def correct(self, tm, h, s, seq, oligo, conc0=0):
        return getattr(self, "%s_correction" % self.__mode.name.lower())(
            seq=seq, h=h, s=s, tm=tm, oligo=oligo, conc0=conc0
        )


class Melter(object):
    DEFAULT_DEGREE_TYPE = const.DEGREE_TYPES.KELVIN
    DEFAULT_NN = "DNA:DNA"
    DEFAULT_OLIGO = 0.25e-6
    DEFAULT_RANGE = 10
    DEFAULT_STEP = 0.1

    __oligo = DEFAULT_OLIGO
    ions = MeltingIonCorrector()
    denaturant = MeltingDenaturantCorrector()
    __nnet = NN_TABLES[DEFAULT_NN]
    __degrees = DEFAULT_DEGREE_TYPE
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
        assert d in const.DEGREE_TYPES
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
        """Calculate melting temperature of provided sequence.
        Args:
            seq (Sequence/string)
        Returns:
            tuple: sequence name, hybridization enthalpy, enthropy, melting
                   temperature, and sequence
        """
        if not type(seq) == Sequence:
            seq = Sequence(seq, const.NATYPES[self.nnet.natypes[0]])
        seq.assert_ab()

        # 1 M NaCl case; SantaLucia, PNAS(95), 1998
        name, g, h, s, tmStd, text = self.__calculate_standard(seq, True)

        # Adjust for FA; Wright, Appl. env. microbiol.(80), 2014
        # Or McConaughy, Biochemistry(8), 1969
        tm = self.denaturant.correct(tmStd, h, s, seq, self.oligo)

        # Adjust for [Na]; Owczarzy et al, Biochemistry(43), 2004
        # Adjust for Mg; Owczarzy et al, Biochemistry(47), 2008
        tm = self.ions.correct(tm, seq)

        g = h - (37 + 273.15) * s
        if self.degrees == const.DEGREE_TYPES.CELSIUS:
            tm -= 273.15
        return (name, g, h, s, tm, text)

    def __calculate_standard(self, seq, forceKelvin=False):
        """Calculate melting temperature of a duplex at standard 1 M NaCl
        (monovalent ions conc). Based on SantaLucia, PNAS(95), 1998
        Args:
            seq (Sequence/string)
        Returns:
            tuple: sequence name, hybridization enthalpy, entropy, melting
                   temperature, and sequence
        """

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

        tm = h / (s + const.R * np.log(self.__oligo))

        g = h - (37 + 273.15) * s
        if self.degrees == const.DEGREE_TYPES.CELSIUS and not forceKelvin:
            tm -= 273.15
        return (seq.name, g, h, s, tm, seq.text)

    def __compute_curve_step_with_mvalue(self, t, h, s, seq):
        m = self.denaturant.mvalue()
        k = self.__dissoc_fraction(t, h, s, m * self.denaturant.conc)
        return (t, k)

    def __compute_curve_step_with_correction(self, t, h, s, seq):
        k = self.__dissoc_fraction(t, h, s, 0)
        t = self.denaturant.correct(t, h, s, seq, self.oligo)
        return (t, k)

    def __run_melting_curve(self, seq, h=None, s=None, tm=None, correctIons=True):
        tmparser = (
            (lambda x: x - 273.15)
            if self.degrees == self.DEGREE_TYPES.CELSIUS
            else (lambda x: x)
        )
        if self.denaturant.mode == self.denaturant.MODES.WRIGHT:
            tm = self.denaturant.correct(tm, h, s, seq, self.oligo)
            curvec = self.__compute_curve_step_with_mvalue
        else:
            curvec = self.__compute_curve_step_with_correction
        curve = []
        tstart = tm - self.curve_range / 2
        if correctIons:
            for ti in range(int(self.curve_range / self.curve_step)):
                t, f = curvec(tstart + self.curve_step * ti, h, s, seq)
                t = self.ions.correct(t, seq)
                curve.append((tmparser(t), f))
        else:
            for ti in range(int(self.curve_range / self.curve_step)):
                t, f = curvec(tstart + self.curve_step * ti, h, s, seq)
                curve.append((tmparser(t), f))
        return curve

    def melting_curve(self, seq, h=None, s=None, tm=None, correctIons=True):
        """Generate melting curve.
        Args:
            seq (Sequence)
            g (float): hybridization free energy
            h (float): hybridization enthalpy
            s (float): hybridization entropy
            tm (float): melting temperature
        Returns:
            list of tuples: [(temperature, dissociated fraction), ...]
        """

        if not type(seq) == Sequence:
            seq = Sequence(seq, const.NATYPES[self.nnet.natypes[0]])
        seq.assert_ab()

        ntype = type(None)
        if ntype == type(h) or ntype == type(s) or ntype == type(tm):
            name, g, h, s, tm, text = self.__calculate_standard(seq, True)

        curve = np.array(self.__run_melting_curve(seq, h, s, tm, correctIons))

        return curve

    def __dissoc_fraction(self, t, h, s, gplus):
        """Calculate duplex dissociated fraction at given temperature.
        Args:
            t (float): temperature
            h (float): enthalpy
            s (float): entropy
            g (float): free energy difference due to denaturants
        Returns:
            float: dissociated fraction
        """
        dg = h - t * s + gplus
        factor = np.exp(-dg / (const.R * t)) * self.oligo
        return 1 / (1 + factor)
