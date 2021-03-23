"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from oligo_melting.asserts import enable_rich_assert
from oligo_melting.const import DEGREE_TYPES, DENATURANT_MODES, NATYPES
from oligo_melting.entab import NN_TABLES
from oligo_melting.melt import Melter, MeltingIonCorrector, MeltingDenaturantCorrector
from oligo_melting.sequence import Sequence
from oligo_melting.scripts import arguments as ap
from Bio.SeqIO.FastaIO import SimpleFastaParser  # type: ignore
from joblib import Parallel, delayed  # type: ignore
import logging
import multiprocessing as mp
import os
import pandas as pd  # type: ignore
from tqdm import tqdm  # type: ignore


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "duplex",
        description="""
Calculate melting temperature of a DNA duplex at provided [oligo],
[Na+], [Mg2+]. Either provide an oligo sequence or a FASTA file.
References:
 [1] Freier et al, PNAS(83), 1986;
 [2] Sugimoto et al, Biochemistry(34), 1995.
 [3] Allawi & Santalucia, Biochemistry(36), 1997;
 [4] SantaLucia, PNAS(95), 1998;
 [5] Owczarzy et al, Biochemistry(43), 2004;
 [6] Owczarzy et al, Biochemistry(47), 2008;
 [7] McConaughy et al, Biochemistry(8), 1969;
 [8] Wright et al, Appl. env. microbiol.(80), 2014.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Calculate melting temperature of NA:NA duplexes.",
    )

    parser.add_argument(
        "input",
        type=str,
        help="""DNA duplex sequence
        (one strand only) or path to a FASTA file.""",
    )
    parser.add_argument(
        "-T",
        type=str,
        help="""Duplex type. Possible values: DNA:DNA (based on ref.3, default),
        RNA:RNA (based on ref.1), DNA:RNA (based on ref.2., given DNA sequence)
        or RNA:DNA (based on ref.2, given RNA sequence). The first nucleic acid type
        indicates the provided sequence.""",
        choices=NN_TABLES.keys(),
        default=Melter.DEFAULT_NN,
    )

    output = parser.add_argument_group("output arguments")
    output.add_argument(
        "--curve-range",
        type=float,
        metavar="range",
        help=f"""Range size to scan around the melting temperature, to build the
        melting curve. Default: {Melter.DEFAULT_RANGE}""",
        default=Melter.DEFAULT_RANGE,
    )
    output.add_argument(
        "--curve-step",
        type=float,
        metavar="range",
        help=f"""Temperature step used to scan the specified range, when building
        the melting curve. Default: {Melter.DEFAULT_STEP}""",
        default=Melter.DEFAULT_STEP,
    )
    output.add_argument(
        "-O", metavar="output", type=str, help="""Path to output file."""
    )

    chemistry = parser.add_argument_group("chemistry-related arguments")
    chemistry.add_argument(
        "-n",
        metavar="na_conc",
        type=float,
        help=f"""Na+ concentration [M].
        Default: {MeltingIonCorrector.DEFAULT_MONOVALENT:.2E} M""",
        default=MeltingIonCorrector.DEFAULT_MONOVALENT,
    )
    chemistry.add_argument(
        "-m",
        metavar="mg_conc",
        type=float,
        help=f"""Mg2+ concentration [M]. Note: Mg2+ correction overwrites Na+
        correction. Default: {MeltingIonCorrector.DEFAULT_DIVALENT:.2E} M""",
        default=MeltingIonCorrector.DEFAULT_DIVALENT,
    )
    chemistry.add_argument(
        "-f",
        type=float,
        metavar="fa_conc",
        help=f"""Formamide concentration in perc.v,v.
        Default: {MeltingDenaturantCorrector.DEFAULT_CONC:.2f}""",
        default=MeltingDenaturantCorrector.DEFAULT_CONC,
    )
    chemistry.add_argument(
        "-o",
        metavar="oligo_conc",
        type=float,
        help=f"""Oligonucleotide concentration [M].
        Default: {Melter.DEFAULT_OLIGO:.2E} M""",
        default=Melter.DEFAULT_OLIGO,
    )
    chemistry.add_argument(
        "--f-mode",
        metavar="fa_mode",
        help="""Mode of formamide correction. "mcconaughy"
        for classical -0.72%%FA correction from ref. 7, "wright" for single reaction
        model correction from ref.8 (default).""",
        choices=DENATURANT_MODES,
        default=MeltingDenaturantCorrector.DEFAULT_MODE,
    )
    chemistry.add_argument(
        "--f-mvalue-1",
        type=str,
        metavar="m1",
        help=f"""Specify the formamide m-value to be used with the wright
        correction model. If the second m-value is 0, this is used as a single value
        "x", otherwise the two values are used as "xL+y" where L is the probe
        length. Default: {MeltingDenaturantCorrector.DEFAULT_M1}""",
        default=MeltingDenaturantCorrector.DEFAULT_M1,
    )
    chemistry.add_argument(
        "--f-mvalue-2",
        type=str,
        metavar="m2",
        help=f"""Specify the formamide m-value to be used with the wright
        correction model. This is used as y in "xL+y", where L is the probe length.
        Default: {MeltingDenaturantCorrector.DEFAULT_M2}""",
        default=MeltingDenaturantCorrector.DEFAULT_M2,
    )

    advanced = parser.add_argument_group("advanced arguments")
    advanced.add_argument(
        "-d",
        "--delim",
        type=str,
        metavar="sep",
        help='''Delimiter between key and value in FASTA header description.
        Default: "="''',
        default="=",
    )
    advanced.add_argument(
        "-t",
        type=int,
        metavar="threads",
        help="""Number of threads for parallelization. Default: 1""",
        default=1,
    )
    parser.add_argument(
        "--out-curve",
        type=str,
        metavar="output",
        help="""Provide path to file where melting curve points should be
        stored.""",
    )
    advanced.add_argument(
        "-C",
        "--celsius",
        dest="celsius",
        action="store_const",
        const=True,
        default=False,
        help="Output temperature in Celsius degrees. Default: Kelvin",
    )
    advanced.add_argument(
        "-F",
        "--fasta-like",
        dest="fasta_like",
        action="store_const",
        const=True,
        default=False,
        help="Output in FASTA format.",
    )

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


def parse_output(args: argparse.Namespace, x):
    if args.fasta_like:
        return ">%s tm%s%.2f\n%s" % (x[0], args.delim, x[4], x[-1])
    else:
        return "%s\t%.3f\t%.3f\t%.3f\t%.2f\t%s" % x


def parse_print_and_write(args: argparse.Namespace, x):
    x = parse_output(args, x)
    print(x)
    args.OH.write("%s\n" % x)


def parse_and_print(args: argparse.Namespace, x):
    x = parse_output(args, x)
    print(x)


def calc_melting_and_curve(record, NATYPE, melter):
    seq = Sequence(record[1], NATYPE, record[0])
    curve = pd.DataFrame(melter.melting_curve(seq))
    curve["name"] = record[0]
    return {"melt": melter.calculate(seq), "curve": curve}


def calc_melting(record, NATYPE, melter):
    seq = Sequence(record[1], NATYPE, record[0])
    return {"melt": melter.calculate(seq)}


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    assert args.o >= 0, "concentration cannot be negative."
    assert args.n >= 0, "concentration cannot be negative."
    assert args.m >= 0, "concentration cannot be negative."
    assert args.f >= 0, "concentration cannot be negative."
    assert args.curve_range > 0, "temperature range must be positive."
    assert args.curve_step > 0, "temperature step must be positive."
    args.t = min(mp.cpu_count(), args.t)

    args.OH = None
    if args.O is not None:
        args.OH = open(args.O, "w+")
        args.pout = parse_print_and_write
    else:
        args.pout = parse_and_print

    if args.out_curve is not None:
        if os.path.isfile(args.out_curve):
            os.remove(args.out_curve)
        args.COH = open(args.out_curve, "a+")

    if args.out_curve is not None:
        args.meltc = calc_melting_and_curve
    else:
        args.meltc = calc_melting

    return args


def init_melter(args: argparse.Namespace) -> Melter:
    args.header = "name\tdG\tdH\tdS\tTm\tSeq"
    melter = Melter()
    melter.oligo = args.o
    melter.ions.monovalent = args.n
    melter.ions.divalent = args.m
    melter.denaturant.conc = args.f
    melter.denaturant.mode = args.f_mode
    melter.denaturant.m1 = args.f_mvalue_1
    melter.denaturant.m2 = args.f_mvalue_2
    melter.load_nn_table(args.T)
    if args.celsius:
        melter.degrees = DEGREE_TYPES.CELSIUS
    melter.curve_range = args.curve_range
    melter.curve_step = args.curve_step
    return melter


def run_single_thread(args, NATYPE, melter):
    if not args.fasta_like:
        print(args.header)
    for record in SimpleFastaParser(args.IH):
        seq = Sequence(record[1], NATYPE, record[0])
        data = melter.calculate(seq)
        args.pout(args, data)
        if args.out_curve is not None:
            curve = pd.DataFrame(melter.melting_curve(seq))
            curve["name"] = data[0]
            curve.to_csv(args.COH, "\t", index=None, header=None)


def run_parallel(args, NATYPE, melter):
    data = Parallel(n_jobs=args.t, verbose=11)(
        delayed(args.meltc)(record, NATYPE, melter)
        for record in SimpleFastaParser(args.IH)
    )
    if not args.fasta_like:
        print(args.header)
    if args.out_curve is not None:
        for record in tqdm(data):
            args.pout(args, record["melt"])
    else:
        for record in tqdm(data):
            args.pout(args, record["melt"])
            record["curve"].to_csv(args.COH, "\t", index=False, header=False)
    return data


def run_single_sequence(args, melter):
    data = melter.calculate(args.input)
    args.pout(args, data)
    if args.out_curve is not None:
        curve = pd.DataFrame(melter.melting_curve(args.input))
        curve["name"] = data[0]
        curve.to_csv(args.COH, "\t", index=None, header=None)
    return data


def close_handles(args):
    if args.OH is not None:
        args.OH.close()
    if args.out_curve is not None:
        args.COH.close()


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    melter = init_melter(args)
    if not args.fasta_like and args.OH is not None:
        args.OH.write("%s\n" % args.header)

    NATYPE = NATYPES[args.T.split(":")[0]]

    if os.path.isfile(args.input):
        args.IH = open(args.input, "r")
        if 1 == args.t:
            run_single_thread(args, NATYPE, melter)
        else:
            run_parallel(args, NATYPE, melter)
        args.IH.close()
    else:
        run_single_sequence(args, melter)

    close_handles(args)
