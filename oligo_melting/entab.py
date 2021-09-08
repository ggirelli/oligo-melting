"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from oligo_melting.const import NATYPES, NN_TABLES_PATH
import os
import pandas as pd  # type: ignore
import pkg_resources  # type: ignore
import re


class NNEnergyTable(object):
    def __init__(self, path, natypes):
        super(NNEnergyTable, self).__init__()

        natypes = natypes.upper()
        assert re.match("[A-Z]+:[A-Z]+", natypes) is not None
        assert all(x in dir(NATYPES) for x in natypes.split(":"))
        self.__natypes = natypes.split(":")

        assert os.path.isfile(path), "'%s' file not found." % path
        self.__table = pd.read_csv(path, "\t", index_col=0, names=["dH0", "dS0", "dG0"])

        self.has_init = "init" in self.__table.index
        if self.has_init:
            self.__init = self.__table.loc["init", :].to_dict()
            self.__table.drop("init", axis=0)

        self.has_sym = "sym" in self.__table.index
        if self.has_sym:
            self.__sym = self.__table.loc["sym", :].to_dict()
            self.__table.drop("sym", axis=0)

        self.has_end = any(x.startswith("end") for x in self.__table.index)
        if self.has_end:
            endRows = [x for x in self.__table.index if x.startswith("end")]
            self.__end = self.__table.loc[endRows, :]
            self.__table.drop(endRows, axis=0)
            self.__end.index = [x[-1] for x in self.__end.index]
            self.__end = self.__end.to_dict()

        self.__ab = set("".join(x for x in self.__table.index if len(x) == 2))
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


PKGROOT = pkg_resources.resource_filename("oligo_melting", "")
NN_TABLES = dict(
    [
        (":".join(x.natypes), x)
        for x in [
            NNEnergyTable(os.path.join(PKGROOT, v), k)
            for k, v in NN_TABLES_PATH.items()
        ]
    ]
)
