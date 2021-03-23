"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from oligo_melting.scripts import arguments
from oligo_melting.scripts import melt, melt_duplex, melt_secstr

import logging
from rich.logging import RichHandler  # type: ignore

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)

__all__ = ["arguments", "melt", "melt_duplex", "melt_secstr"]
