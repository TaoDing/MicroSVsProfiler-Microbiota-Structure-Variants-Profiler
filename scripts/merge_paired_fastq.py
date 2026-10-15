#!/usr/bin/env python3
"""Convert paired FASTQ files to a single FASTQ stream for HUMAnN.

Pairs are checked for matching identifiers and emitted consecutively (R1 then R2).
This does not assemble or merge overlapping paired-end sequences.
"""
import argparse
import gzip
from pathlib import Path


def reader(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as f:
        while True:
            record = [f.readline() for _ in range(4)]
            if not record[0]:
                if any(record[1:]):
                    raise ValueError(f'Incomplete trailing FASTQ record: {path}')
                return
            if not all(record) or not record[0].startswith('@') or not record[2].startswith('+'):
                raise ValueError(f'Invalid FASTQ record: {path}')
            if len(record[1].rstrip('\r\n')) != len(record[3].rstrip('\r\n')):
                raise ValueError(f'Sequence/quality lengths differ: {path}')
            yield record


def identifier(header):
    name = header[1:].split()[0]
    return name.removesuffix('/1').removesuffix('/2')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--r1', required=True)
    p.add_argument('--r2', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    from itertools import zip_longest
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    try:
        with output.open('w', encoding='utf-8') as out:
            for mate1, mate2 in zip_longest(reader(args.r1), reader(args.r2)):
                if mate1 is None or mate2 is None:
                    raise ValueError('Paired FASTQ files contain different numbers of records')
                if identifier(mate1[0]) != identifier(mate2[0]):
                    raise ValueError(f'Unmatched pair at record {n+1}: {mate1[0].strip()}, {mate2[0].strip()}')
                out.writelines(mate1)
                out.writelines(mate2)
                n += 1
    except Exception:
        output.unlink(missing_ok=True)
        raise
    print(f'Validated {n} read pairs; wrote {2*n} reads to {output}')


if __name__ == '__main__':
    main()
