#!/usr/bin/env python3
"""Convert trusted SGV-Finder pandas pickles to CSV (not for untrusted pickles)."""
import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dsv', required=True, type=Path, help='Trusted dsgv.df pickle')
    parser.add_argument('--vsv', required=True, type=Path, help='Trusted vsgv.df pickle')
    parser.add_argument('--out', default='results/step3/converted', type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for source, target in ((args.dsv, 'dsgv.csv'), (args.vsv, 'vsgv.csv')):
        if not source.is_file():
            parser.error(f'File not found: {source}')
        df = pd.read_pickle(source)  # Only load files you created/trust.
        df.to_csv(args.out / target)
        print(f'Saved {args.out / target} ({df.shape[0]} rows x {df.shape[1]} columns)')


if __name__ == '__main__':
    main()
