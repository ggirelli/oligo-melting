"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from oligo_melting import nearest_neighbor

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    pass

__all__ = ["nearest_neighbor"]
