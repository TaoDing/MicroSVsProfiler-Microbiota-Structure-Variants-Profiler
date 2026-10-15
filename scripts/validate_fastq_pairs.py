#!/usr/bin/env python3
"""Stream paired FASTQ files; validate four-line records and mate IDs."""
import gzip
import sys
from itertools import zip_longest


def records(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, plus, q = fh.readline(), fh.readline(), fh.readline()
            if not s or not plus or not q or not h.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'Malformed FASTQ record in {path}')
            if len(s.rstrip('\r\n')) != len(q.rstrip('\r\n')):
                raise ValueError(f'Length mismatch in {path}: {h.strip()}')
            yield h[1:].split()[0].removesuffix('/1').removesuffix('/2')


def main():
    if len(sys.argv) != 3:
        sys.exit('Usage: validate_fastq_pairs.py R1.fastq[.gz] R2.fastq[.gz]')
    count = 0
    for count, (id1, id2) in enumerate(zip_longest(records(sys.argv[1]), records(sys.argv[2])), 1):
        if id1 is None or id2 is None or id1 != id2:
            sys.exit(f'Paired FASTQ mismatch at record {count}: {id1} vs {id2}')
    print(f'Validated {count} read pairs')


if __name__ == '__main__':
    main()
