"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from oligo_melting import __version__
import sys

# import multiprocessing as mp


def add_version_option(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--version", action="version", version=f"{sys.argv[0]} {__version__}"
    )
    return parser
