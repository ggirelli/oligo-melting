"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from oligo_melting import __version__
from oligo_melting.scripts import arguments as ap
from oligo_melting import scripts
import sys


def default_parser(*args) -> None:
    print("melt -h for usage details.")
    sys.exit()


def main():
    parser = argparse.ArgumentParser(
        description=f"""
Version:    {__version__}
Author:     Gabriele Girelli
Docs:       http://ggirelli.github.io/oligo_melting
Code:       http://github.com/ggirelli/oligo_melting

A Python3 package for melting temperature calculation
of oligonucleotides hybridization and secondary structures.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.set_defaults(parse=default_parser)
    parser = ap.add_version_option(parser)

    subparsers = parser.add_subparsers(
        title="sub-commands",
        help="Access the help page for a sub-command with: sub-command -h",
    )

    scripts.melt_duplex.init_parser(subparsers)
    scripts.melt_secstr.init_parser(subparsers)

    args = parser.parse_args()
    args = args.parse(args)
    args.run(args)
