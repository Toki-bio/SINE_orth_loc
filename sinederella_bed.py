#!/usr/bin/env python3
"""SINE copies of a SINEderella run as the BED input of SINE_orth_loc.bash in bank mode.

Usage: sinederella_bed.py RUN_DIR -o sp-BANK.bed [--no-partial]

Takes from RUN_DIR (a SINEderella result directory):
  loci_verdict.tsv(.gz)              copies with verdict SINE (assigned, flanks unique); grey and not_SINE are left out
  family_view/partial_copies.tsv     partial copies (Status "partial"): label <best consensus>:<a>-<b>, the span of
                                     the best consensus they cover (their verdict is not_SINE only because they have
                                     no subfamily; their flanks are tested by the flank mapping of SINE_orth_loc)
Copy names are chrom:start-end(strand) with BED starts and '_' written as '@U@' (SINEderella genome cleaning).
Overlapping copies on one strand become one locus with the label of the longest copy (a whole composite, not
its pieces). Writes chrom, start, end, label, 0, strand, sorted; counts to stderr. Python 3 standard library."""
import argparse, gzip, os, re, sys

NAME = re.compile(r'^(.*):(\d+)-(\d+)\(([+-])\)$')


def opener(path):
    if os.path.exists(path):
        return open(path, encoding='utf-8')
    if os.path.exists(path + '.gz'):
        return gzip.open(path + '.gz', 'rt', encoding='utf-8')
    sys.exit(f'sinederella_bed.py: {path}(.gz) not found')


def parse(name):
    m = NAME.match(name.replace('@U@', '_'))
    if not m:
        sys.exit(f'sinederella_bed.py: cannot read copy name {name}')
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run')
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--no-partial', action='store_true', help='leave partial copies out')
    a = ap.parse_args()

    rows, n = [], {'SINE': 0, 'partial': 0}
    with opener(os.path.join(a.run, 'loci_verdict.tsv')) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        ic, isf, iv = head.index('copy'), head.index('subfamily'), head.index('verdict')
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[iv] == 'SINE':
                rows.append(parse(f[ic]) + (f[isf], 'SINE'))
                n['SINE'] += 1
    if not a.no_partial:
        with opener(os.path.join(a.run, 'family_view', 'partial_copies.tsv')) as fh:
            head = fh.readline().rstrip('\n').split('\t')
            ic, ib, ist, ista = head.index('Sequence'), head.index('Best'), head.index('Stretch'), head.index('Status')
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if f[ista] != 'partial':
                    continue
                span = [s for s in f[ist].split(',') if s.rsplit(':', 1)[0] == f[ib]]
                if not span:
                    sys.exit(f'sinederella_bed.py: no span of {f[ib]} in {f[ist]}')
                rows.append(parse(f[ic]) + (span[0], 'partial'))
                n['partial'] += 1

    # one locus per group of overlapping copies on one strand: the longest copy's label and coordinates span all
    rows.sort(key=lambda r: (r[0], r[3], r[1], r[2]))
    out, merged = [], 0
    for r in rows:
        if out and out[-1][0] == r[0] and out[-1][3] == r[3] and r[1] < out[-1][2]:
            g = out[-1]
            keep = r if (r[2] - r[1]) > (g[6] - g[5]) else None
            out[-1] = [g[0], g[1], max(g[2], r[2]), g[3], keep[4] if keep else g[4], keep[1] if keep else g[5],
                       keep[2] if keep else g[6]]
            merged += 1
        else:
            out.append([r[0], r[1], r[2], r[3], r[4], r[1], r[2]])
    out.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(a.out, 'w', encoding='utf-8', newline='\n') as fo:
        for c, s, e, st, lab, _, _ in out:
            fo.write(f'{c}\t{s}\t{e}\t{lab}\t0\t{st}\n')
    print(f'{a.run}: {n["SINE"]} SINE + {n["partial"]} partial copies; {merged} overlapping copies merged; '
          f'{len(out)} loci -> {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
