"""Command line entry point for the primary analysis.

Reads one raw chain from ``data/`` and writes the result tables into
``results/``. Run it from the repository root::

    python python_scripts/analyze_and_save.py --bhw 10 --nt 200 \
        --nstep 1000000 --skip 50000

The numerics live in :mod:`analysis.pipeline`; this file only maps command
line arguments onto the directory layout.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.pipeline import analyze_dataset  # noqa: E402


def build_parser():
    """Return the argument parser for this script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nt', type=int, required=True,
                        help='number of time slices of the run to analyse')
    parser.add_argument('--bhw', type=int, required=True,
                        help='beta * hbar * omega of the run')
    parser.add_argument('--nstep', type=int, default=1000000,
                        help='number of measurements in the raw file')
    parser.add_argument('--skip', type=int, default=0,
                        help='leading measurements to drop as thermalisation')
    return parser


def main(argv=None):
    """Analyse the dataset selected by the command line arguments."""
    args = build_parser().parse_args(argv)

    basedir = f'bhw{args.bhw}_nstep{args.nstep}'
    raw_path = os.path.join('data', basedir, f'raw_data_nt{args.nt}.dat')
    outdir = os.path.join('results', basedir,
                          f'nt{args.nt}_therm{args.skip}')
    eta = float(args.bhw) / args.nt

    print(f'nt={args.nt}, eta={eta:.6f}')
    analyze_dataset(raw_path, outdir, eta, skip=args.skip, verbose=True)


if __name__ == '__main__':
    main()
