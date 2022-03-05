"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from enum import Enum

R = 1.987 / 1000  # kcal / (K mol)


class DEGREE_TYPES(Enum):
    CELSIUS = 1
    KELVIN = 2


class DENATURANT_MODES(Enum):
    MCCONAUGHY = 1
    WRIGHT = 2


AB_DNA = ["ACGT", "TGCA"]
AB_RNA = ["ACGU", "UGCA"]


class NATYPES(Enum):
    DNA = 1
    RNA = 2


AB_NA = {NATYPES.DNA: AB_DNA, NATYPES.RNA: AB_RNA}


NN_TABLES_PATH = {
    # Table from Freier et al, PNAS(83), 1986 - in 1 M NaCl [RNA]
    "RNA:RNA": "../nntables/freier.tsv",
    # Table from Sugimoto et al, Biochemistry(34), 1995 - in 1 M NaCl [DNA:RNA]
    # For calculation from the DNA sequence 5'-to-3'
    "DNA:RNA": "../nntables/sugimoto_2.tsv",
    # Table from Sugimoto et al, Biochemistry(34), 1995 - in 1 M NaCl [DNA:RNA]
    # For calculation from the RNA sequence 5'-to-3'
    "RNA:DNA": "../nntables/sugimoto_1.tsv",
    # Table from Allawi&Santalucia, Biochemistry(36), 1997 - in 1 M NaCl [DNA]
    "DNA:DNA": "../nntables/allawy.tsv",
}
