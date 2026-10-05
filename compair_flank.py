#!/usr/bin/env python3
"""Helpers for ComPair.sh (Python standard library only).

  compair_flank.py mask ALN FROM TO [MINRUN]
      Identity of the first two sequences of alignment ALN over columns FROM..TO (1-based, inclusive)
      after removing columns that belong to one-sided gap runs of >= MINRUN (default 50) in either
      sequence: a large insertion/deletion next to the insertion site (often another element) is set
      aside instead of failing the flank. Prints "PID NID NCOLS NMASKED" with esl-alipid's convention
      PID = 100 * identities / min(residues of seq1, residues of seq2) over the kept columns.

  compair_flank.py edge ALN SIZES
      Prints the 0-based index of the first of the two loci whose window touches the start (0) or the
      end of its sequence, according to SIZES (a .fai or chrom<TAB>size file of both genomes); prints
      -1 if neither does. A window clipped at a contig end means the sequence is missing, not diverged.
"""
import re
import sys

LOCUS = re.compile(r'^>?(.+):(\d+)-(\d+)\([+-]\)')


def read_aln(path, n=2):
    names, seqs = [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if len(names) == n:
                    break
                names.append(line[1:]); seqs.append([])
            elif seqs:
                seqs[-1].append(line)
    return names, [''.join(s) for s in seqs]


def mask(path, start, end, minrun=50):
    _, (a, b) = read_aln(path)
    a, b = a[start - 1:end], b[start - 1:end]
    drop = [False] * len(a)
    for s in (a, b):
        for m in re.finditer('-{%d,}' % minrun, s):
            for i in range(m.start(), m.end()):
                drop[i] = True
    ra = rb = ident = cols = 0
    for x, y, d in zip(a, b, drop):
        if d:
            continue
        cols += 1
        ra += x != '-'; rb += y != '-'
        ident += x != '-' and x.upper() == y.upper() and x.upper() != 'N'
    pid = 100.0 * ident / min(ra, rb) if min(ra, rb) else 0.0
    return pid, ident, cols, sum(drop)


def edge(path, sizes_path):
    sizes = {}
    with open(sizes_path) as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    names, _ = read_aln(path)
    for i, name in enumerate(names):
        m = LOCUS.match(name.split('::')[0])
        if not m:
            continue
        chrom, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        if s <= 0 or (chrom in sizes and e >= sizes[chrom]):
            return i
    return -1


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    if sys.argv[1] == 'mask':
        minrun = int(sys.argv[5]) if len(sys.argv) > 5 else 50
        pid, nid, cols, masked = mask(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), minrun)
        print(f'{int(pid)} {nid} {cols} {masked}')
    elif sys.argv[1] == 'edge':
        print(edge(sys.argv[2], sys.argv[3]))
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
