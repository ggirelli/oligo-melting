"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from oligo_melting.const import AB_NA


class Sequence(object):
    __len = None
    __fgc = None
    __rc = None

    def __init__(self, seq, t, name=None):
        super(Sequence, self).__init__()
        self.__text = seq.upper()
        self.__len = len(self.__text)
        self.__natype = t
        self.__ab = AB_NA[t]
        self.__name = "%d-mer" % self.__len if name is None else name

    @property
    def text(self):
        return self.__text

    @property
    def len(self):
        return self.__len

    @property
    def rc(self):
        if self.__rc is None:
            self.__rc = self.mkrc(self.text, self.natype)
        return self.__rc

    @property
    def fgc(self):
        if self.__fgc is None:
            if self.__len != 0:
                self.__fgc = self.__text.count("G")
                self.__fgc += self.__text.count("C")
                self.__fgc /= self.__len
            else:
                self.__fgc = 0
        return self.__fgc

    @property
    def natype(self):
        return self.__natype

    @property
    def ab(self):
        return self.__ab

    @property
    def name(self):
        return self.__name

    def __eq__(self, other):
        if any(attrname not in dir(other) for attrname in ["text", "name", "natype"]):
            return False
        if not self.text == other.text:
            return False
        if not self.name == other.name:
            return False
        if not self.natype == other.natype:
            return False
        return True

    def dimers(self):
        """Extract NN dimers from sequence.
        Args:
            seq (str): nucleic acid sequence.
        Returns:
            list: NN dimers.
        """
        return self.dimerator(self.text)

    def assert_ab(self):
        assert_msg = "sequence alphabet and nucleic acid type mismatch."
        assert_msg += "\n%s\t%s" % (set(self.text), self.ab[0])
        assert self.check_ab(self.text, self.ab), assert_msg

    @staticmethod
    def check_ab(seq, ab):
        return all(x in ab[0] for x in set(seq))

    @staticmethod
    def mkrc(na, t):
        """Calculate reverse complement.
        Args:
            na (string): nucleic acid sequence.
            t (string): nucleic acid type, either 'dna' or 'rna'.
        Return:
            string: reverse complement of na.
        """

        assert t in AB_NA.keys(), "%s not in %s" % (t, list(AB_NA.keys()))
        ab = AB_NA[t]

        rab = ab[1].strip().lower()
        ab = ab[0].strip().lower()
        na = na.lower()

        for c in na:
            assert_msg = "provided string conflicts with selected alphabet."
            assert_msg += " (%s in %s)" % (c, na)
            assert c in ab, assert_msg

        r = na[::-1]
        rc = [rab[ab.index(c)] for c in r]
        rc = "".join(str(c) for c in rc).upper()

        return rc

    @staticmethod
    def dimerator(seq):
        """Extract NN dimers from sequence.
        Args:
            seq (str): nucleic acid sequence.
        Returns:
            list: NN dimers.
        """
        return (seq[i : (i + 2)] for i in range(len(seq) - 1))
