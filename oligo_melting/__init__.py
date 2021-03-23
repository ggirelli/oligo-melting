"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from oligo_melting import asserts, const
from oligo_melting import entab, melt, sequence

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    pass

__all__ = ["asserts", "const", "entab", "melt", "sequence"]
